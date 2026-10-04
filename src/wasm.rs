//! C-ABI exports for the browser build (`wasm32-unknown-unknown`, no bindings generator):
//! js/engine.js copies program text into memory from `ft_alloc`, builds a program, and reads
//! sampled records back from the returned buffer.

use crate::batch::{sample, Backend};
use crate::noise::{FaultPlan, SiteTree};
use crate::program::Program;
use std::cell::RefCell;

pub struct WasmProgram {
    program: Program,
    sites: SiteTree,
    out: Vec<u8>,
    dist: Vec<f64>,
}

thread_local! {
    static LAST_ERROR: RefCell<Vec<u8>> = const { RefCell::new(Vec::new()) };
}

fn set_error(e: String) {
    LAST_ERROR.with(|l| *l.borrow_mut() = e.into_bytes());
}

#[no_mangle]
pub extern "C" fn ft_alloc(len: usize) -> *mut u8 {
    let mut v = Vec::<u8>::with_capacity(len.max(1));
    let p = v.as_mut_ptr();
    std::mem::forget(v);
    p
}

/// # Safety
/// `ptr` must come from `ft_alloc(len)` and not be freed twice.
#[no_mangle]
pub unsafe extern "C" fn ft_free(ptr: *mut u8, len: usize) {
    drop(Vec::from_raw_parts(ptr, 0, len.max(1)));
}

#[no_mangle]
pub extern "C" fn ft_error_ptr() -> *const u8 {
    LAST_ERROR.with(|l| l.borrow().as_ptr())
}

#[no_mangle]
pub extern "C" fn ft_error_len() -> usize {
    LAST_ERROR.with(|l| l.borrow().len())
}

/// Parses program text (UTF-8, `len` bytes at `ptr`); null on error (see `ft_error_*`).
///
/// # Safety
/// `ptr..ptr+len` must be readable.
#[no_mangle]
pub unsafe extern "C" fn ft_program_new(ptr: *const u8, len: usize) -> *mut WasmProgram {
    let bytes = std::slice::from_raw_parts(ptr, len);
    let text = match std::str::from_utf8(bytes) {
        Ok(t) => t,
        Err(_) => {
            set_error("program text is not UTF-8".into());
            return std::ptr::null_mut();
        }
    };
    match Program::parse(text) {
        Ok(program) => {
            let sites = SiteTree::build(&program);
            Box::into_raw(Box::new(WasmProgram { program, sites, out: Vec::new(), dist: Vec::new() }))
        }
        Err(e) => {
            set_error(e);
            std::ptr::null_mut()
        }
    }
}

/// # Safety
/// `p` must come from `ft_program_new` and not be freed twice.
#[no_mangle]
pub unsafe extern "C" fn ft_program_free(p: *mut WasmProgram) {
    if !p.is_null() {
        drop(Box::from_raw(p));
    }
}

/// # Safety
/// `p` must be a live program.
#[no_mangle]
pub unsafe extern "C" fn ft_num_measurements(p: *const WasmProgram) -> u32 {
    (*p).program.num_measurements as u32
}

/// # Safety
/// `p` must be a live program.
#[no_mangle]
pub unsafe extern "C" fn ft_expected_faults(p: *const WasmProgram) -> f64 {
    (*p).sites.expected_faults()
}

/// [P(K=0) … P(K=kmax), P(K>kmax)]; returns a pointer to kmax + 2 f64s.
///
/// # Safety
/// `p` must be a live program.
#[no_mangle]
pub unsafe extern "C" fn ft_fault_distribution(p: *mut WasmProgram, kmax: u32) -> *const f64 {
    let w = &mut *p;
    w.dist = w.sites.fault_count_distribution(kmax as usize);
    w.dist.as_ptr()
}

/// Samples shots first_shot … first_shot + shots − 1 with stream seed. `faults`: −1 plain,
/// −2 none, k ≥ 0 exactly k. Returns a pointer to shots × ⌈m/8⌉ bytes (null on error).
///
/// # Safety
/// `p` must be a live program.
#[no_mangle]
pub unsafe extern "C" fn ft_sample(p: *mut WasmProgram, shots: u32, seed: f64, first_shot: f64, faults: i32) -> *const u8 {
    let w = &mut *p;
    let plan = match faults {
        -1 => FaultPlan::Plain,
        -2 => FaultPlan::None,
        k if k >= 0 => FaultPlan::Exactly(k as usize),
        _ => {
            set_error("bad fault mode".into());
            return std::ptr::null();
        }
    };
    match sample(&w.program, &w.sites, shots as u64, seed as u64, first_shot as u64, plan, Backend::Sparse, 1, false) {
        Ok(r) => {
            w.out = r.records;
            w.out.as_ptr()
        }
        Err(e) => {
            set_error(e);
            std::ptr::null()
        }
    }
}

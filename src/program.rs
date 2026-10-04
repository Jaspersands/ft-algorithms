//! Logical programs: the text format, its parser and static validation.
//!
//! One instruction per line; `#` starts a comment. Blocks are `REPEAT n { … }` and
//! `TABLE rec[-a] rec[-b] … { case 0 } { case 1 } …` (2^k cases chosen by the named records,
//! the first named record being bit 0 of the case index).

use crate::gates::{Basis, Gate, Noise};

pub type Block = Vec<Op>;

#[derive(Clone, Debug, PartialEq)]
pub enum Op {
    Gate(Gate, Vec<u32>),
    Measure(Basis, Vec<u32>),
    Reset(Basis, Vec<u32>),
    Noise(Noise, Vec<u32>),
    Mark(u16),
    Repeat(u64, Block),
    /// Record offsets (k ≥ 1 means rec[-k]) and the 2^offsets.len() cases.
    Table(Vec<u32>, Vec<Block>),
}

#[derive(Clone, Debug)]
pub struct Program {
    pub body: Block,
    pub num_qubits: u32,
    pub num_measurements: u64,
    /// Fault classes named by MARK; index 0 is "unmarked".
    pub classes: Vec<String>,
}

/// Limits that keep a hostile program from exhausting memory or time before it runs.
pub const MAX_QUBIT: u32 = 127;
pub const MAX_MEASUREMENTS: u64 = 1 << 32;
pub const MAX_TABLE_BITS: usize = 16;

struct Line {
    no: usize,
    text: String,
}

impl Program {
    pub fn parse(text: &str) -> Result<Program, String> {
        // Split braces onto their own lines, keeping the source line number of each piece.
        let mut lines = Vec::new();
        for (i, raw) in text.lines().enumerate() {
            let no = i + 1;
            let code = match raw.find('#') {
                Some(k) => &raw[..k],
                None => raw,
            };
            let mut cur = String::new();
            for ch in code.chars() {
                if ch == '{' || ch == '}' {
                    if !cur.trim().is_empty() {
                        lines.push(Line { no, text: cur.trim().to_string() });
                    }
                    cur.clear();
                    lines.push(Line { no, text: ch.to_string() });
                } else {
                    cur.push(ch);
                }
            }
            if !cur.trim().is_empty() {
                lines.push(Line { no, text: cur.trim().to_string() });
            }
        }
        let mut parser = Parser { lines, pos: 0, classes: vec!["unmarked".to_string()], max_qubit: None };
        let body = parser.block(false)?;
        let mut records = 0u64;
        check_block(&body, &mut records)?;
        let num_qubits = parser.max_qubit.map(|q| q + 1).unwrap_or(0);
        Ok(Program { body, num_qubits, num_measurements: records, classes: parser.classes })
    }
}

struct Parser {
    lines: Vec<Line>,
    pos: usize,
    classes: Vec<String>,
    max_qubit: Option<u32>,
}

fn err(no: usize, msg: impl std::fmt::Display) -> String {
    format!("line {no}: {msg}")
}

impl Parser {
    fn block(&mut self, nested: bool) -> Result<Block, String> {
        let mut ops = Vec::new();
        while self.pos < self.lines.len() {
            let no = self.lines[self.pos].no;
            let text = self.lines[self.pos].text.clone();
            self.pos += 1;
            if text == "}" {
                if nested {
                    return Ok(ops);
                }
                return Err(err(no, "unmatched '}'"));
            }
            if text == "{" {
                return Err(err(no, "'{' without REPEAT or TABLE"));
            }
            self.instruction(no, &text, &mut ops)?;
        }
        if nested {
            return Err("unterminated block: missing '}'".to_string());
        }
        Ok(ops)
    }

    fn open_brace(&mut self, no: usize) -> Result<(), String> {
        if self.pos < self.lines.len() && self.lines[self.pos].text == "{" {
            self.pos += 1;
            Ok(())
        } else {
            Err(err(no, "expected '{'"))
        }
    }

    fn qubits(&mut self, no: usize, toks: &[&str]) -> Result<Vec<u32>, String> {
        let mut out = Vec::with_capacity(toks.len());
        for t in toks {
            let q: u32 = t.parse().map_err(|_| err(no, format!("bad qubit target '{t}'")))?;
            if q > MAX_QUBIT {
                return Err(err(no, format!("qubit {q} exceeds the limit {MAX_QUBIT}")));
            }
            self.max_qubit = Some(self.max_qubit.map_or(q, |m| m.max(q)));
            out.push(q);
        }
        Ok(out)
    }

    fn instruction(&mut self, no: usize, text: &str, ops: &mut Block) -> Result<(), String> {
        // NAME(args) targets…
        let (head, rest) = match text.find(|c: char| c.is_whitespace()) {
            Some(k) if !text[..k].contains('(') || text[..k].contains(')') => (&text[..k], text[k..].trim()),
            _ => match text.find(')') {
                Some(k) => (&text[..=k], text[k + 1..].trim()),
                None => (text, ""),
            },
        };
        let (name, args): (&str, Option<&str>) = match head.find('(') {
            Some(k) => {
                if !head.ends_with(')') {
                    return Err(err(no, format!("malformed '{head}'")));
                }
                (&head[..k], Some(&head[k + 1..head.len() - 1]))
            }
            None => (head, None),
        };
        let toks: Vec<&str> = rest.split_whitespace().collect();
        let nums = |s: Option<&str>| -> Result<Vec<f64>, String> {
            match s {
                None => Ok(vec![]),
                Some(s) => s
                    .split(',')
                    .map(|x| x.trim().parse::<f64>().map_err(|_| err(no, format!("bad number '{x}'"))))
                    .collect(),
            }
        };
        match name {
            "TICK" => Ok(()),
            "MARK" => {
                let label = args.ok_or_else(|| err(no, "MARK needs a name: MARK(name)"))?.trim().to_string();
                if label.is_empty() {
                    return Err(err(no, "empty MARK name"));
                }
                let id = match self.classes.iter().position(|c| *c == label) {
                    Some(i) => i,
                    None => {
                        self.classes.push(label);
                        self.classes.len() - 1
                    }
                };
                if id > u16::MAX as usize {
                    return Err(err(no, "too many MARK classes"));
                }
                ops.push(Op::Mark(id as u16));
                Ok(())
            }
            "REPEAT" => {
                if toks.len() != 1 {
                    return Err(err(no, "REPEAT takes one count"));
                }
                let n: u64 = toks[0].parse().map_err(|_| err(no, format!("bad REPEAT count '{}'", toks[0])))?;
                if n == 0 {
                    return Err(err(no, "REPEAT count must be at least 1"));
                }
                self.open_brace(no)?;
                let body = self.block(true)?;
                ops.push(Op::Repeat(n, body));
                Ok(())
            }
            "TABLE" => {
                if toks.is_empty() || toks.len() > MAX_TABLE_BITS {
                    return Err(err(no, format!("TABLE takes 1 to {MAX_TABLE_BITS} records")));
                }
                let mut offs = Vec::new();
                for t in &toks {
                    let k = t
                        .strip_prefix("rec[-")
                        .and_then(|s| s.strip_suffix(']'))
                        .and_then(|s| s.parse::<u32>().ok())
                        .filter(|&k| k >= 1)
                        .ok_or_else(|| err(no, format!("bad record target '{t}' (want rec[-k], k ≥ 1)")))?;
                    offs.push(k);
                }
                let want = 1usize << offs.len();
                let mut cases = Vec::with_capacity(want);
                for _ in 0..want {
                    self.open_brace(no).map_err(|_| err(no, format!("TABLE with {} records needs {want} case blocks", offs.len())))?;
                    cases.push(self.block(true)?);
                }
                if self.pos < self.lines.len() && self.lines[self.pos].text == "{" {
                    return Err(err(no, format!("TABLE with {} records needs exactly {want} case blocks", offs.len())));
                }
                ops.push(Op::Table(offs, cases));
                Ok(())
            }
            _ => {
                if let Some(g) = Gate::from_name(name) {
                    if args.is_some() {
                        return Err(err(no, format!("{name} takes no arguments")));
                    }
                    let qs = self.qubits(no, &toks)?;
                    let a = g.arity();
                    if qs.is_empty() || qs.len() % a != 0 {
                        return Err(err(no, format!("{name} needs a multiple of {a} targets")));
                    }
                    if a > 1 {
                        for grp in qs.chunks(a) {
                            for i in 0..a {
                                for j in i + 1..a {
                                    if grp[i] == grp[j] {
                                        return Err(err(no, format!("{name} applied to qubit {} twice", grp[i])));
                                    }
                                }
                            }
                        }
                    }
                    ops.push(Op::Gate(g, qs));
                    return Ok(());
                }
                let basis = match name {
                    "M" | "R" => Some(Basis::Z),
                    "MX" | "RX" => Some(Basis::X),
                    "MY" | "RY" => Some(Basis::Y),
                    _ => None,
                };
                if let Some(b) = basis {
                    if args.is_some() {
                        return Err(err(no, format!("{name} takes no arguments")));
                    }
                    let qs = self.qubits(no, &toks)?;
                    if qs.is_empty() {
                        return Err(err(no, format!("{name} needs targets")));
                    }
                    if name.starts_with('M') {
                        ops.push(Op::Measure(b, qs));
                    } else {
                        ops.push(Op::Reset(b, qs));
                    }
                    return Ok(());
                }
                let a = nums(args)?;
                let noise = match (name, a.len()) {
                    ("PAULI_CHANNEL_1", 3) => Noise::P1([a[0], a[1], a[2]]),
                    ("PAULI_CHANNEL_2", 15) => {
                        let mut p = [0.0; 15];
                        p.copy_from_slice(&a);
                        Noise::P2(p)
                    }
                    ("Z_CHANNEL_3", 7) => {
                        let mut p = [0.0; 7];
                        p.copy_from_slice(&a);
                        Noise::Z3(p)
                    }
                    ("X_ERROR", 1) => Noise::P1([a[0], 0.0, 0.0]),
                    ("Y_ERROR", 1) => Noise::P1([0.0, a[0], 0.0]),
                    ("Z_ERROR", 1) => Noise::P1([0.0, 0.0, a[0]]),
                    ("DEPOLARIZE1", 1) => Noise::P1([a[0] / 3.0; 3]),
                    ("DEPOLARIZE2", 1) => Noise::P2([a[0] / 15.0; 15]),
                    ("PAULI_CHANNEL_1" | "PAULI_CHANNEL_2" | "Z_CHANNEL_3" | "X_ERROR" | "Y_ERROR" | "Z_ERROR"
                    | "DEPOLARIZE1" | "DEPOLARIZE2", n) => {
                        return Err(err(no, format!("{name} got {n} arguments")));
                    }
                    _ => return Err(err(no, format!("unknown instruction '{name}'"))),
                };
                for &p in noise.probs() {
                    if !(0.0..=1.0).contains(&p) || p.is_nan() {
                        return Err(err(no, format!("probability {p} outside [0, 1]")));
                    }
                }
                if noise.total() > 1.0 + 1e-12 {
                    return Err(err(no, format!("probabilities sum to {} > 1", noise.total())));
                }
                if noise.total() >= 1.0 {
                    return Err(err(no, "a channel that always fires (total probability 1) is not supported"));
                }
                let qs = self.qubits(no, &toks)?;
                let ar = noise.arity();
                if qs.is_empty() || qs.len() % ar != 0 {
                    return Err(err(no, format!("{name} needs a multiple of {ar} targets")));
                }
                if ar > 1 {
                    for grp in qs.chunks(ar) {
                        for i in 0..ar {
                            for j in i + 1..ar {
                                if grp[i] == grp[j] {
                                    return Err(err(no, format!("{name} applied to qubit {} twice", grp[i])));
                                }
                            }
                        }
                    }
                }
                ops.push(Op::Noise(noise, qs));
                Ok(())
            }
        }
    }
}

/// Number of measurement records one execution of `block` appends.
pub fn block_measurements(block: &Block) -> u64 {
    let mut n = 0u64;
    for op in block {
        n = n.saturating_add(match op {
            Op::Measure(_, qs) => qs.len() as u64,
            Op::Repeat(k, b) => k.saturating_mul(block_measurements(b)),
            Op::Table(_, cases) => block_measurements(&cases[0]),
            _ => 0,
        });
    }
    n
}

/// Checks record references and TABLE case consistency; `records` is the count before `block`.
fn check_block(block: &Block, records: &mut u64) -> Result<(), String> {
    for op in block {
        match op {
            Op::Measure(_, qs) => *records += qs.len() as u64,
            Op::Repeat(k, b) => {
                let start = *records;
                let mut r = start;
                check_block(b, &mut r)?;
                let per = r - start;
                *records = start.saturating_add(per.saturating_mul(*k));
            }
            Op::Table(offs, cases) => {
                for &k in offs {
                    if k as u64 > *records {
                        return Err(format!("TABLE reads rec[-{k}] but only {} records exist there", *records));
                    }
                }
                let m0 = block_measurements(&cases[0]);
                for c in cases {
                    if block_measurements(c) != m0 {
                        return Err("TABLE cases must contain the same number of measurements".to_string());
                    }
                }
                for c in cases {
                    let mut r = *records;
                    check_block(c, &mut r)?;
                }
                *records += m0;
            }
            _ => {}
        }
        if *records > MAX_MEASUREMENTS {
            return Err(format!("more than {MAX_MEASUREMENTS} measurement records"));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_every_form() {
        let p = Program::parse(
            "# comment\nH 0 1\nCX 0 1 2 3\nCCZ 0 1 2\nT_DAG 4\nM 0 1\nMX 2\nR 3\nRX 3\n\
             PAULI_CHANNEL_1(0.1, 0, 0.2) 0\nDEPOLARIZE2(0.01) 0 1\nZ_CHANNEL_3(0.001,0,0,0,0,0,0.002) 0 1 2\n\
             MARK(cnot)\nREPEAT 3 {\n  X 0\n  M 0\n}\nTABLE rec[-1] rec[-2] {X 2} {Y 2} {\nZ 2\n} {H 2}\nTICK\n",
        )
        .unwrap();
        assert_eq!(p.num_qubits, 5);
        assert_eq!(p.num_measurements, 2 + 1 + 3);
        assert_eq!(p.classes, vec!["unmarked".to_string(), "cnot".to_string()]);
        match &p.body[p.body.len() - 1] {
            Op::Table(offs, cases) => {
                assert_eq!(offs, &vec![1, 2]);
                assert_eq!(cases.len(), 4);
                assert_eq!(cases[2], vec![Op::Gate(Gate::Z, vec![2])]);
            }
            other => panic!("{other:?}"),
        }
        match &p.body[8] {
            Op::Noise(Noise::P1(pr), qs) => {
                assert_eq!(pr, &[0.1, 0.0, 0.2]);
                assert_eq!(qs, &vec![0]);
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn nested_repeat_counts() {
        let p = Program::parse("REPEAT 4 {\nREPEAT 5 {\nM 0 1\n}\nM 2\n}").unwrap();
        assert_eq!(p.num_measurements, 4 * (10 + 1));
    }

    fn bad(text: &str, needle: &str) {
        let e = Program::parse(text).unwrap_err();
        assert!(e.contains(needle), "{text:?} gave {e:?}, wanted {needle:?}");
    }

    #[test]
    fn rejects_invalid() {
        bad("FOO 0", "unknown instruction");
        bad("CX 0", "multiple of 2");
        bad("CX 1 1", "twice");
        bad("CCZ 0 1 1", "twice");
        bad("H", "multiple of 1");
        bad("H x", "bad qubit");
        bad("H 128", "exceeds");
        bad("X_ERROR(1.5) 0", "outside");
        bad("PAULI_CHANNEL_1(0.5,0.5,0.1) 0", "sum");
        bad("PAULI_CHANNEL_1(0.1) 0", "got 1 arguments");
        bad("REPEAT 0 {\nX 0\n}", "at least 1");
        bad("REPEAT 2\nX 0", "expected '{'");
        bad("REPEAT 2 {\nX 0", "unterminated");
        bad("}", "unmatched");
        bad("M 0\nTABLE rec[-2] {X 0} {X 1}", "only 1 records");
        bad("M 0\nTABLE rec[-1] {X 0}", "needs 2 case blocks");
        bad("M 0\nTABLE rec[-1] {X 0} {X 1} {X 2}", "exactly 2");
        bad("M 0\nTABLE rec[-1] {M 0} {X 1}", "same number of measurements");
        bad("M 0\nTABLE rec[0] {X 0} {X 1}", "bad record target");
        bad("H(0.1) 0", "takes no arguments");
        bad("MARK 0", "MARK needs a name");
        bad("X_ERROR(0.1) 0\nX_ERROR(1) 0", "always fires");
    }

    #[test]
    fn line_numbers_in_errors() {
        let e = Program::parse("H 0\n\nCX 0 0").unwrap_err();
        assert!(e.starts_with("line 3:"), "{e}");
    }

    #[test]
    fn table_in_repeat_sees_body_records() {
        Program::parse("REPEAT 3 {\nM 0\nTABLE rec[-1] {X 1} {Z 1}\n}").unwrap();
    }
}

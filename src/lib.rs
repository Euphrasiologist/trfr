// Requires trf

use core::fmt;
use std::{
    error::Error as StdError,
    fs::File,
    io::{self, BufRead, Read},
    num::{ParseFloatError, ParseIntError},
    path::Path,
    result::Result as StdResult,
};

pub type Result<T> = StdResult<T, Error>;

#[derive(Debug)]
pub struct Error(Box<ErrorKind>);

impl Error {
    /// A crate private constructor for `Error`.
    pub(crate) fn new(kind: ErrorKind) -> Error {
        Error(Box::new(kind))
    }

    /// Return the specific type of this error.
    pub fn kind(&self) -> &ErrorKind {
        &self.0
    }

    /// Unwrap this error into its underlying type.
    pub fn into_kind(self) -> ErrorKind {
        *self.0
    }
}

#[derive(Debug)]
pub enum ErrorKind {
    Io(io::Error),
    Int(ParseIntError),
    Float(ParseFloatError),
    Parser(String),
}

impl From<io::Error> for Error {
    fn from(err: io::Error) -> Error {
        Error::new(ErrorKind::Io(err))
    }
}

impl From<ParseIntError> for Error {
    fn from(err: ParseIntError) -> Self {
        Error::new(ErrorKind::Int(err))
    }
}
impl From<ParseFloatError> for Error {
    fn from(err: ParseFloatError) -> Self {
        Error::new(ErrorKind::Float(err))
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self.0 {
            ErrorKind::Io(ref err) => write!(f, "I/O error - {}", err),
            ErrorKind::Int(ref err) => write!(f, "parsing integer error - {}", err),
            ErrorKind::Float(ref err) => write!(f, "parsing float error - {}", err),
            ErrorKind::Parser(ref err) => write!(f, "parser error - {}", err),
        }
    }
}

impl StdError for Error {}

#[derive(Debug, Default, PartialEq)]
pub struct Record {
    // the name of the fasta record
    seq_id: String,
    // start index of the repeat
    start: usize,
    // end index of the repeat
    end: usize,
    // the period of the repeat
    period: u16,
    // number of copies aligned with the consensus pattern
    copy_number: f32,
    // size of consensus pattern (may differ slightly from the period size).
    consensus_pattern_size: u16,
    // percent of matches between adjacent copies overall
    perc_matches: u8,
    // percent of indels between adjacent copies overall
    perc_indels: u8,
    // alignment score
    alignment_score: u32,
    // percentages of the nucleotides
    perc_a: u8,
    perc_c: u8,
    perc_g: u8,
    perc_t: u8,
    // entropy measure based on percent composition
    entropy: f32,
    // the repeat pattern itself
    consensus_pattern: String,
    // and the longer repeat sequence extracted from the reads
    repeat_seq: String,
}

pub struct Reader<R> {
    rdr: io::BufReader<R>,
    line: u64,
    id: String,
}

impl Reader<File> {
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Reader<File>> {
        Ok(Reader::new(File::open(path)?))
    }
}

impl<R: io::Read> Reader<R> {
    pub fn new(rdr: R) -> Reader<R> {
        Reader {
            rdr: io::BufReader::new(rdr),
            line: 0,
            id: String::new(),
        }
    }

    /// A borrowed iterator over the records of a refer file.
    pub fn records(&mut self) -> RecordsIter<R> {
        RecordsIter::new(self)
    }

    /// An owned iterator over the records of a refer file.
    pub fn into_records(self) -> RecordsIntoIter<R> {
        RecordsIntoIter::new(self)
    }

    /// Read a single record from an input reader.
    fn read_record(&mut self) -> Result<Option<Record>> {
        // read the line
        // if line is empty/only contains spaces, end read
        // if line contains something,
        let mut record = Record::default();
        let reader = self.rdr.by_ref();

        // might have to read each line separately
        // using read_line()
        let mut temp_buf = String::new();

        loop {
            self.line += 1;
            temp_buf.clear();
            let bytes = reader.read_line(&mut temp_buf)?;
            if bytes == 0 {
                // this is the EOF
                if record == Record::default() {
                    return Ok(None);
                } else {
                    return Ok(Some(record));
                }
            }

            let parsed = parse_input_line(temp_buf.clone(), &mut record);

            match parsed {
                Ok(e) => match e {
                    // we got something
                    Some(sequence_id) => {
                        if let Some(id) = sequence_id {
                            self.id = id;
                        }
                        break;
                    }
                    // we got a text line which wasn't data
                    None => continue,
                },
                Err(e) => {
                    return Err(Error::new(ErrorKind::Parser(format!(
                        "at line {}, {}",
                        self.line, e
                    ))))
                }
            }
        }

        Ok(Some(record))
    }
}

fn parse_input_line(input: String, record: &mut Record) -> Result<Option<Option<String>>> {
    //

    if input.trim().is_empty() {
        return Ok(None);
    }

    if [
        "Tandem Repeats",
        "Gary Benson",
        "Program",
        "Boston",
        "Version",
        "Parameters",
    ]
    .iter()
    .any(|s| input.starts_with(*s))
    {
        return Ok(None);
    }

    if input.starts_with("Sequence") {
        let name = input.replace("Sequence: ", "");
        let name = name.trim().to_string();
        return Ok(Some(Some(name)));
    }

    // that should cover everything? Now we split the line
    let line_elements = input.split(' ').collect::<Vec<&str>>();

    if let [start, end, period, copy_number, consensus_pattern_size, perc_matches, perc_indels, alignment_score, perc_a, perc_c, perc_g, perc_t, entropy, consensus_pattern, repeat_seq] =
        &line_elements[..]
    {
        record.start = start.parse::<usize>()?;
        record.end = end.parse::<usize>()?;
        record.period = period.parse::<u16>()?;
        record.copy_number = copy_number.parse::<f32>()?;
        record.consensus_pattern_size = consensus_pattern_size.parse::<u16>()?;
        record.perc_matches = perc_matches.parse::<u8>()?;
        record.perc_indels = perc_indels.parse::<u8>()?;
        record.alignment_score = alignment_score.parse::<u32>()?;
        record.perc_a = perc_a.parse::<u8>()?;
        record.perc_c = perc_c.parse::<u8>()?;
        record.perc_g = perc_g.parse::<u8>()?;
        record.perc_t = perc_t.parse::<u8>()?;
        record.entropy = entropy.parse::<f32>()?;
        record.consensus_pattern = consensus_pattern.to_string();
        record.repeat_seq = repeat_seq.trim().to_string();

        return Ok(Some(None));
    }

    Ok(None)
}

/// A borrowed iterator over the records of a refer file.
pub struct RecordsIter<'r, R: 'r> {
    /// The underlying reader
    rdr: &'r mut Reader<R>,
}

impl<'r, R: io::Read> RecordsIter<'r, R> {
    fn new(rdr: &'r mut Reader<R>) -> RecordsIter<'r, R> {
        RecordsIter { rdr }
    }
    /// Return a reference to the underlying reader.
    pub fn reader(&self) -> &Reader<R> {
        self.rdr
    }

    /// Return a mutable reference to the underlying reader.
    pub fn reader_mut(&mut self) -> &mut Reader<R> {
        self.rdr
    }
}

impl<'r, R: io::Read> Iterator for RecordsIter<'r, R> {
    type Item = Result<Record>;

    fn next(&mut self) -> Option<Result<Record>> {
        match self.rdr.read_record() {
            Ok(Some(mut r)) => {
                self.rdr.line += 1;
                r.seq_id = self.rdr.id.clone();
                Some(Ok(r))
            }
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

/// An owned iterator over the records of a refer file.
pub struct RecordsIntoIter<R> {
    /// The underlying reader.
    rdr: Reader<R>,
}

impl<R: io::Read> RecordsIntoIter<R> {
    fn new(rdr: Reader<R>) -> RecordsIntoIter<R> {
        RecordsIntoIter { rdr }
    }
    /// Return a reference to the underlying reader.
    pub fn reader(&self) -> &Reader<R> {
        &self.rdr
    }

    /// Return a mutable reference to the underlying reader.
    pub fn reader_mut(&mut self) -> &mut Reader<R> {
        &mut self.rdr
    }

    /// Drop this iterator and return the underlying reader.
    pub fn into_reader(self) -> Reader<R> {
        self.rdr
    }
}

impl<R: io::Read> Iterator for RecordsIntoIter<R> {
    type Item = Result<Record>;

    fn next(&mut self) -> Option<Result<Record>> {
        match self.rdr.read_record() {
            Ok(Some(mut r)) => {
                self.rdr.line += 1;
                r.seq_id = self.rdr.id.clone();
                Some(Ok(r))
            }
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

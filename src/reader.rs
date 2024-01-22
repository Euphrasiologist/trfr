use std::{
    fs::File,
    io::{self, BufRead, Read},
    path::Path,
};

use crate::{Error, ErrorKind, Record, Result};

#[derive(Clone, Copy)]
pub enum Flag {
    D,
    Ngs,
}

pub struct Reader<R> {
    rdr: io::BufReader<R>,
    line: u64,
    id: String,
    flag: Flag,
}

impl Reader<File> {
    pub fn from_path<P: AsRef<Path>>(path: P, flag: Flag) -> Result<Reader<File>> {
        Ok(Reader::new(File::open(path)?, flag))
    }
    pub fn from_reader<R: io::Read>(rdr: R, flag: Flag) -> Reader<R> {
        Reader::new(rdr, flag)
    }
}

impl<R: io::Read> Reader<R> {
    pub fn new(rdr: R, flag: Flag) -> Reader<R> {
        Reader {
            rdr: io::BufReader::new(rdr),
            line: 0,
            id: String::new(),
            flag,
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

            let parsed = parse_input_line(temp_buf.clone(), &mut record, self.flag);

            match parsed {
                Ok(e) => match e {
                    // we got something
                    Some(sequence_id) => {
                        if let Some(id) = sequence_id {
                            self.id = id;
                            continue;
                        }
                        break;
                    }
                    // we got a text line which wasn't data
                    None => continue,
                },
                Err(e) => {
                    return Err(Error::new(ErrorKind::ReadRecord(format!(
                        "at line {}, {}",
                        self.line, e
                    ))))
                }
            }
        }

        Ok(Some(record))
    }
}

/// An inner private function to parse an input line. The output type is a little
/// complex as we want to either skip lines, save the sequence name, or append
/// fields to the `Record` struct.
fn parse_input_line(
    input: String,
    record: &mut Record,
    flag: Flag,
) -> Result<Option<Option<String>>> {
    // depending on the flag
    match flag {
        // aka -d
        Flag::D => {
            // do the stuff we were doing before
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
                let name = input.replace("Sequence: ", "").trim().to_string();
                return Ok(Some(Some(name)));
            }
        }
        // aka -ngs
        Flag::Ngs => {
            // simpler here, as we only worry about @ lines
            if input.starts_with('@') {
                let name = input.replace('@', "").trim().to_string();
                return Ok(Some(Some(name)));
            }
        }
    }

    // that should cover everything? Now we split the line
    let mut line_elements = input.split(' ').collect::<Vec<&str>>();
    // ignore flanking regions
    match flag {
        Flag::D => (),
        Flag::Ngs => {
            // two extra columns, so remove them
            line_elements.truncate(line_elements.len() - 2);
            assert!(line_elements.len() == 15);
        }
    }

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

        Ok(Some(None))
    } else {
        Err(Error::new(ErrorKind::Parser(
            "could not split into 15 elements".into(),
        )))
    }
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

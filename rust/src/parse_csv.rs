use csv::{Reader, StringRecord};
use na::{Vector3, Quaternion, Unit};
use serde::Deserialize;
use serde_with::with_prefix;
use std::{error::Error, fs};
use nalgebra as na;
use itertools::Itertools;
use std::slice;

pub type Result<T> = std::result::Result<T, Box<dyn Error>>;

with_prefix!(prefix_patella "patella ");
with_prefix!(prefix_probe "probe ");
with_prefix!(prefix_pin1 "pin1 ");
with_prefix!(prefix_pin2 "pin2 ");

#[derive(Deserialize, Debug)]
#[allow(dead_code)]
pub struct Datum {
    frame: u32,
    #[serde(flatten, with = "prefix_probe")]
    probe: Option<DatumCoords>,
    #[serde(flatten, with = "prefix_patella")]
    patella: Option<DatumCoords>,
    #[serde(flatten, with = "prefix_pin1")]
    pin1: Option<DatumCoords>,
    #[serde(flatten, with = "prefix_pin2")]
    pin2: Option<DatumCoords>,
}

#[derive(Deserialize, Clone, Copy, Debug)]
pub struct DatumCoords {
#[serde(deserialize_with = "csv::invalid_option")]
    q0: Option<f64>,
#[serde(deserialize_with = "csv::invalid_option")]
    qx: Option<f64>,
#[serde(deserialize_with = "csv::invalid_option")]
    qy: Option<f64>,
#[serde(deserialize_with = "csv::invalid_option")]
    qz: Option<f64>,
#[serde(deserialize_with = "csv::invalid_option")]
    x: Option<f64>,
#[serde(deserialize_with = "csv::invalid_option")]
    y: Option<f64>,
#[serde(deserialize_with = "csv::invalid_option")]
    z: Option<f64>,
#[serde(deserialize_with = "csv::invalid_option")]
    error: Option<f64>,
}

#[derive(Debug)]
pub struct Trackers {
    pub probe: Record,
    pub pin1: Record,
    pub pin2: Record,
}

#[allow(dead_code)]
#[derive(Debug, Default)]
pub struct Record {
    pub q0: Vec<f64>,
    pub qx: Vec<f64>,
    pub qy: Vec<f64>,
    pub qz: Vec<f64>,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub z: Vec<f64>,
    pub error: Vec<f64>,
}

pub trait Mean {
    fn mean(&self) -> f64;
}
impl Mean for Vec<f64> {
    fn mean(&self) -> f64 {
        self.iter().sum::<f64>()/(self.len() as f64)
    }
}

#[derive(Deserialize, Debug)]
#[serde(rename_all = "lowercase")]
pub enum Units {
    Mm,
    Cm,
}

pub struct Csv {
    pub header: Header,
    pub content: Trackers,
}

#[allow(dead_code)]
#[derive(Deserialize, Debug)]
#[serde(rename_all = "PascalCase")]
pub struct Header {
    #[serde(rename = "Number of frames")]
    pub number_of_frames: usize,
    pub frequency: f64,
    pub units: Units,
}

impl Csv {
    pub fn new(name: &str) -> Result<Self> {
        let raw = fs::read_to_string(name)?;
        let mut headers = raw.lines().take(5);
        let header: Header = serde_yaml::from_str(&headers.by_ref().take(3).join("\n"))?;
        // let header: Header = serde_yaml::from_str(&raw.lines().take(3).collect::<Vec<&str>>().join("\n"))?;

        let csv_headers: Vec<String> = headers
            .by_ref()
            .last()
            .unwrap()
            .split(',')
            .map(|header| match header.rsplit_once(char::is_whitespace) {
                Some((name, coord)) if !name.contains("Pin") => format!("probe {}", coord.to_lowercase()),
                _ => header.to_lowercase(),
            })
        .collect();
        
        
        let csv_content: String = raw.lines()
            .skip(5)
            .map(|f| f.chars().filter(|c| !c.is_whitespace()).collect::<String>())
            .join("\n");

        let mut rdr = Reader::from_reader(csv_content.as_bytes());
        rdr.set_headers(StringRecord::from(csv_headers));

        let mut content = vec![];
        for result in rdr.deserialize() {
            content.push(result?);
        }
        let content: Trackers = Self::convert_to_structure(content);

        Ok(Csv { header, content })
    }

    fn convert_to_structure(vec_datum: Vec<Datum>) -> Trackers {
        let temp = vec_datum.iter();
        let probe = Record::new(&temp, |f| f.probe);
        let pin1 = Record::new(&temp, |f| f.pin1);
        let pin2 = Record::new(&temp, |f| f.pin2);

        Trackers { probe, pin1, pin2 }
    }
}

// #[allow(dead_code)]
impl Record {
    pub fn new<F>(data: &slice::Iter<Datum>, structure: F) -> Record
        where
            F: Fn(&Datum) -> Option<DatumCoords>,
    {
        let (q0, qx, qy, qz, x, y, z, error): (Vec<_>,Vec<_>,Vec<_>,Vec<_>,Vec<_>,Vec<_>,Vec<_>,Vec<_>) =
           data.as_ref().iter()
           .filter_map(structure)
           .map(|f| (f.q0.unwrap_or(0.0), f.qx.unwrap_or(0.0), f.qy.unwrap_or(0.0), f.qz.unwrap_or(0.0), f.x.unwrap_or(0.0), f.y.unwrap_or(0.0), f.z.unwrap_or(0.0), f.error.unwrap_or(0.0)))
           .multiunzip();
        Record { q0, qx, qy, qz, x, y, z, error}
    }
    // pub fn mean(&self) -> Coords {
    //     Coords {
    //         q0: vec![self.q0.mean()],
    //         qx: vec![self.qx.mean()],
    //         qy: vec![self.qy.mean()],
    //         qz: vec![self.qz.mean()],
    //         x: vec![self.x.mean()],
    //         y: vec![self.y.mean()],
    //         z: vec![self.z.mean()],
    //         error: vec![self.error.mean()]
    //     }
    // }

    pub fn rotation(c: (&f64,&f64,&f64,&f64)) -> Unit<Quaternion<f64>> {
        Unit::new_normalize(Quaternion::new(*c.0, *c.1, *c.2, *c.3))
    }

    // pub fn tracker<R: Reference>(&self) -> Vec<CoordinateSystem<R>> {
    //     izip!(&self.q0, &self.qx, &self.qy, &self.qz)
    //         .map(Self::rotation)
    //         .map(Matrix3::from)
    //         .enumerate()
    //         .map(|(i, rotation)| {
    //             let mut transform = na::Matrix4::identity();
    //
    //             transform.fixed_view_mut::<3,3>(0, 0).copy_from(&rotation);
    //             transform.fixed_view_mut::<3,1>(0, 3).copy_from(&Vector3::new(self.x[i], self.y[i], self.z[i]));
    //             CoordinateSystem::from_transform(transform)
    //         })
    //         .collect()
    // }

    pub fn mean_xyz(&self) -> Vector3<f64> {
        Vector3::new(self.x.mean(), self.y.mean(), self.z.mean())
    }

}

use std::marker::PhantomData;

use crate::{anatomy::RigidBodyLandmarks, parse_csv::Coords};

use super::anatomy::{Bone, Landmark, Position, Side};
use itertools::izip;
use nalgebra::{
    self as na, Isometry3, Matrix3, Matrix4, Rotation3, Translation3, Vector3,
};
// use super::change_frame::ChangeFrame;

pub trait Reference: Sized {}

pub trait Inverse<R: Reference, T: Reference> {
    fn inverse(&self) -> Vec<CoordinateSystem<T, R>>
    where
        Self: Sized;
}

// pub struct Global<'a> {
//     pub pin1: &'a [Matrix4<f64>],
//     pub pin2: &'a [Matrix4<f64>],
// }
#[derive(Debug)]
pub struct Global;
#[derive(Debug)]
pub struct Tibia;
#[derive(Debug)]
pub struct Femur;

#[allow(dead_code)]
#[derive(Debug)]
pub struct Patella;
#[derive(Debug)]
pub struct Pin1;
#[derive(Debug)]
pub struct Pin2;

#[allow(dead_code)]
pub struct UnitVectors {
    pub x: Vector3<f64>,
    pub y: Vector3<f64>,
    pub z: Vector3<f64>,
}

#[derive(Debug)]
pub struct CoordinateSystem<O, G>
where
    O: Reference,
    G: Reference,
{
    pub transform: Isometry3<f64>,
    reference: PhantomData<O>,
    in_reference_to: PhantomData<G>,
}

#[derive(Debug)]
struct Tracker<O: Reference, G: Reference> {
    coordinate_system: CoordinateSystem<O, G>,
}

// impl<'a, O: Reference, G: Reference> Tracker<O, G> {
//     pub fn new<N: Reference> (coords: &Coords) -> Tracker<N, Global>{
//         Tracker { coordinate_system: () }
//     }
// }

impl<O, G> CoordinateSystem<O, G>
where
    O: Reference,
    G: Reference,
{
    pub fn from_isometry<N: Reference>(transform: Isometry3<f64>) -> CoordinateSystem<N, G> {
        CoordinateSystem {
            transform,
            reference: PhantomData,
            in_reference_to: PhantomData,
        }
    }
    pub fn from_transform<N: Reference>(transform: Matrix4<f64>) -> CoordinateSystem<N, G> {
        let translation = Vector3::from(transform.fixed_view::<3, 1>(0, 3));
        let rotation = Rotation3::from_matrix(&transform.fixed_view::<3, 3>(0, 0).into());
        let isometry = Isometry3::from_parts(translation.into(), rotation.into());
        CoordinateSystem {
            transform: isometry,
            reference: PhantomData,
            in_reference_to: PhantomData,
        }
    }
    pub fn change_frame<N>(
        &self,
        matrix: &CoordinateSystem<G, N>,
    ) -> CoordinateSystem<O, N>
    where
        N: Reference,
    {
        self.transform(matrix)
    }
    pub fn origin(&self) -> &Translation3<f64> {
        &self.transform.translation
    }
    pub fn tracker<N: Reference>(coords: &Coords) -> Vec<CoordinateSystem<N, G>> {
        izip!(&coords.q0, &coords.qx, &coords.qy, &coords.qz)
            .map(Coords::rotation)
            .map(Matrix3::from)
            .enumerate()
            .map(|(i, rotation)| {
                let mut transform = Matrix4::identity();

                transform.fixed_view_mut::<3, 3>(0, 0).copy_from(&rotation);
                transform
                    .fixed_view_mut::<3, 1>(0, 3)
                    .copy_from(&Vector3::new(coords.x[i], coords.y[i], coords.z[i]));
                CoordinateSystem::<O, G>::from_transform::<N>(transform)
            })
            .collect()
    }
    pub fn floating_axis<N: Reference>(&self, other: &CoordinateSystem<N, G>) -> Vector3<f64> {
        self.unit_vectors()
            .z
            .cross(&other.unit_vectors().x)
            .normalize()
    }
    // pub fn from_landmark(anatomy: &[Landmark]) -> CoordinateSystem<O, Global> {
    //     let (transform, point) = body_frame(anatomy);
    //     CoordinateSystem {
    //         transform,
    //         origin: Vector4::new(point[0], point[1], point[2], 1.0),
    //         reference: PhantomData::<O>,
    //         in_reference_to: PhantomData::<Global>,
    //     }
    // }
    pub fn rotation(&self, side: &Side) -> Vector3<f64> {
        let ang  = self.transform.rotation.euler_angles();
        match side {
            Side::Left => Vector3::new(-ang.0, -ang.1, -ang.2),
            Side::Right => Vector3::new(-ang.0, ang.1, ang.2),
        }
    }
    pub fn unit_vectors(&self) -> UnitVectors {
        let f = self.transform.to_homogeneous();
        UnitVectors {
            x: f.fixed_view::<3, 1>(0, 0).into_owned(),
            y: f.fixed_view::<3, 1>(0, 1).into_owned(),
            z: f.fixed_view::<3, 1>(0, 2).into_owned(),
        }
    }
    // pub fn apply<X: Reference>(
    //     &self,
    //     matrix: &CoordinateSystem<G, X>,
    //     offset: &Matrix4<f64>,
    // ) -> CoordinateSystem<O, X> {
    //     CoordinateSystem {
    //         transform: matrix.transform * self.transform * offset,
    //         reference: PhantomData,
    //         in_reference_to: PhantomData,
    //     }
    // }
    pub fn inverse(&self) -> CoordinateSystem<G, O> {
        CoordinateSystem::<G, O>::from_isometry(self.transform.inverse())
    }
    pub fn transform<X: Reference>(
        &self,
        matrix: &CoordinateSystem<G, X>,
    ) -> CoordinateSystem<O, X> {
        CoordinateSystem {
            transform: matrix.transform * self.transform,
            reference: PhantomData,
            in_reference_to: PhantomData,
        }
    }
}

// impl <O: Reference, G: Reference> FrameChanger<O,G> for Option<CoordinateSystem<O, G>> {
//     fn change_frame<X: Reference>(&self) -> Option<ChangeFrame<O, X ,G>> {
//         match self {
//             Some(old) => Some(ChangeFrame::<O,X,G>::new(old)),
//             _ => None,
//         }
//     }
// }

impl<O, G> Inverse<O, G> for Vec<CoordinateSystem<O, G>>
where
    O: Reference,
    O: Reference,
    G: Reference,
{
    fn inverse(&self) -> Vec<CoordinateSystem<G, O>>
    where
        Self: Sized,
    {
        let invert = |f: &CoordinateSystem<O, G>| {
            CoordinateSystem::<G, O>::from_isometry(f.transform.inverse())
        };
        self.iter().map(invert).collect()
    }
}

pub fn coordinate_system_from_landmark<O: Reference>(
    anatomy: &[Landmark],
) -> CoordinateSystem<O, Global> {
    let transform = body_frame(anatomy);
    CoordinateSystem {
        transform,
        reference: PhantomData::<O>,
        in_reference_to: PhantomData::<Global>,
    }
}

fn body_frame(bone_locations: &[Landmark]) -> Isometry3<f64> {
    let mut locations = RigidBodyLandmarks::default();

    for location in bone_locations.iter() {
        match location.position[0] {
            Position::Medial => locations.medial = Some(location),
            Position::Lateral => locations.lateral = Some(location),
            Position::Anterior => locations.anterior = Some(location),
            Position::Posterior => locations.posterior = Some(location),
            Position::Distal => locations.distal = Some(location),
            Position::Proximal => locations.proximal = Some(location),
        }
    }
    let med = locations
        .medial
        .map(|f| &f.data.probe)
        .map(|f| f.mean_xyz())
        .unwrap();
    let lat = locations
        .lateral
        .map(|f| &f.data.probe)
        .map(|f| f.mean_xyz())
        .unwrap();

    let origin = (med + lat) / 2.0;
    let i = match bone_locations[0].side {
        Side::Right => (lat - med).normalize(),
        Side::Left => (med - lat).normalize(),
    };
    // The distal (tibia) and proximal (femur) points are approximate, thus this axis is not
    // necessarily perpendicular to epicondylar axis.
    let tempk = match bone_locations[0].bone {
        Bone::Tibia => {
            origin
                - (locations
                    .distal
                    .map(|f| &f.data.probe)
                    .map(|f| f.mean_xyz())
                    .unwrap())
        }
        Bone::Femur => {
            locations
                .proximal
                .map(|f| &f.data.probe)
                .map(|f| f.mean_xyz())
                .unwrap()
                - origin
        }
        Bone::Patella => todo!(),
    };

    let j = (tempk.normalize()).cross(&i);
    let k = i.cross(&j).normalize();
    let rotation = Rotation3::from_matrix(&Matrix3::from_columns(&[i, j, k]));
    Isometry3::from_parts(origin.into(), rotation.into())

    // let mut matrix = Matrix4::identity();
    // // rot.fixed_view_mut::<3, 3>(0, 0)
    // //     .copy_from(&Matrix3::from_columns(&[i, j, k]));
    // matrix
    //     .fixed_view_mut::<3, 4>(0, 0)
    //     .copy_from(&Matrix3x4::from_columns(&[i, j, k, origin]));
    //
    // matrix
}
// impl <'a> Reference for Global <'a> { }
impl Reference for Global {}
impl Reference for Tibia {}
impl Reference for Femur {}
impl Reference for Patella {}
impl Reference for Pin1 {}
impl Reference for Pin2 {}

use std::f64::consts::PI;

use nalgebra::Vector3;
use serde::Serialize;

use crate::{
    anatomy::Side,
    frame_of_reference::{CoordinateSystem, Reference},
};

#[derive(Serialize, Clone, Debug)]
pub struct Output {
    pub flexion: f64,
    pub tibial_varus: f64,
    pub tibial_external_rotation: f64,
    pub lateral_medial: f64,
    pub anterior_posterior: f64,
    pub distal_proximal: f64,
}
impl Output {
    pub fn new<O: Reference, N: Reference, G: Reference>(
        first: &CoordinateSystem<O, G>,
        other: &CoordinateSystem<N, G>,
        side: &Side,
    ) -> Output {
        let self_to_other: CoordinateSystem<O, N> = first
            // .change_frame::<Femur>()
            .transform::<N>(&other.inverse());
        let rotation = self_to_other.rotation(side);

        let self_vecs = first.unit_vectors();
        let other_vec = other.unit_vectors();

        let floating_axis = first.floating_axis(other);
        let femur_to_tibia_origin = (first.origin().to_homogeneous()
            - other.origin().to_homogeneous())
        .fixed_view::<3, 1>(0, 0)
        .into_owned();

        let (varus, external_rotation, lat_med) = {
            let beta = other_vec.x.dot(&self_vecs.z).acos();
            match side {
                Side::Right => (
                    (-floating_axis).dot(&self_vecs.x).asin(),
                    (PI / 2.0 - beta),
                    femur_to_tibia_origin.dot(&other_vec.x),
                ),
                Side::Left => (
                    (floating_axis).dot(&self_vecs.x).asin(),
                    -(PI / 2.0 - beta),
                    femur_to_tibia_origin.dot(&-other_vec.x),
                ),
            }
        };
        let ant_post = femur_to_tibia_origin.dot(&floating_axis);
        let dist_prox = -femur_to_tibia_origin.dot(&self_vecs.z);
        let tf_ang = Vector3::new(rotation.x, varus, external_rotation);
        let tf_transl = Vector3::new(lat_med, ant_post, dist_prox);
        let output = Output {
            flexion: tf_ang.x.to_degrees(),
            tibial_varus: tf_ang.y.to_degrees(),
            tibial_external_rotation: tf_ang.z.to_degrees(),
            lateral_medial: tf_transl.x,
            anterior_posterior: tf_transl.y,
            distal_proximal: tf_transl.z,
        };
        println!("{:#?}", output);
        output
    }
}

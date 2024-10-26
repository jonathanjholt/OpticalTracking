use std::f64::consts::PI;

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
        let femur_to_tibia_origin = first.origin().vector - other.origin().vector;

        let (varus, external_rotation, lateral_medial) = {
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
        let anterior_posterior = femur_to_tibia_origin.dot(&floating_axis);
        let distal_proximal = -femur_to_tibia_origin.dot(&self_vecs.z);
        Output {
            flexion: rotation.x.to_degrees(),
            tibial_varus: varus.to_degrees(),
            tibial_external_rotation: external_rotation.to_degrees(),
            lateral_medial,
            anterior_posterior,
            distal_proximal,
        }
    }
}

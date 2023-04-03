//This code has not yet been tested. It may not be correct.

mod gyrovector {
    use nalgebra as na;

    type Point = na::Vector3<f64>;

    pub fn gyrovector_add(u: Point, v: Point, k: f64) -> Point {
        let uv = u.dot(&v);
    
        if k.abs() < f64::EPSILON {
            return u + v;
        }
    
        let denominator = 1.0 - k * (2.0 * uv + uv * uv);
        u + (1.0 / k + uv) * ((1.0 / k + uv + uv * uv) / denominator) * (v - u)
    }
    
    // Gyrovector scalar multiplication with curvature constant K
    pub fn gyrovector_scalar_mul(scalar: f64, v: Point, k: f64) -> Point {
        if k.abs() < f64::EPSILON {
            return scalar * v;
        }
    
        let norm_v = v.norm();
        let lambda = (scalar * norm_v).tanh();
        (v / norm_v) * lambda
    }
    
    // Gyrovector scalar division with curvature constant K
    pub fn gyrovector_scalar_div(v: Point, scalar: f64, k: f64) -> Point {
        if k.abs() < f64::EPSILON {
            return v / scalar;
        }
    
        let norm_v = v.norm();
        let lambda = (norm_v / scalar).tanh();
        (v / norm_v) * lambda
    }
    
    // Gyrovector inversion with curvature constant K
    pub fn gyrovector_inv(v: Point, k: f64) -> Point {
        if k.abs() < f64::EPSILON {
            return -1.0 * v;
        }
    
        gyrovector_scalar_mul(-1.0, v, k)
    }
    
    // Gyrovector dot product with curvature constant K
    pub fn gyrovector_dot(a: Point, b: Point, k: f64) -> f64 {
        let x = a.x * b.x - k * a.y * b.y;
        let y = a.x * b.y + a.y * b.x;
    
        let na = (1.0 - k * a.x * a.x - k * a.y * a.y).sqrt();
        let nb = (1.0 - k * b.x * b.x - k * b.y * b.y).sqrt();
    
        x / (na * nb) - y
    }
    
    // Gyrovector norm with curvature constant K
    pub fn gyrovector_norm(v: Point, k: f64) -> f64 {
        (1.0 - k * v.x * v.x - k * v.y * v.y).sqrt()
    }
    
    // Gyrovector distance with curvature constant K
    pub fn gyrovector_distance(a: Point, b: Point, k: f64) -> f64 {
        let diff = gyrovector_add(a, gyrovector_inv(b, k), k);
        gyrovector_norm(diff, k)
    }

    // Gyrovector cross product with curvature constant K
    pub fn gyrovector_cross(a: Point, b: Point, k: f64) -> Point {
        let euclidean_cross = a.cross(&b);

        if k.abs() < f64::EPSILON {
            return euclidean_cross;
        }

        let norm_a = gyrovector_norm(a, k);
        let norm_b = gyrovector_norm(b, k);

        let coeff = 1.0 / (1.0 - k * (a.dot(&b) / (norm_a * norm_b)).powi(2));
        coeff * euclidean_cross
    }

    // Gyrovector linear interpolation with curvature constant K
    pub fn gyrovector_lerp(a: Point, b: Point, t: f64, k: f64) -> Point {
        let d = gyrovector_distance(a, b, k);
        let s = t * d;
        let direction = gyrovector_scalar_mul(1.0 / d, gyrovector_add(b, gyrovector_inv(a, k), k), k);
        gyrovector_add(a, gyrovector_scalar_mul(s, direction, k), k)
    }

    // Gyrovector rotation with curvature constant K
    pub fn gyrovector_rotate(point: Point, axis: Point, angle: f64, k: f64) -> Point {
        let normalized_axis = axis.normalize();
        let rot_vector = gyrovector_scalar_mul(angle.sinh(), normalized_axis, k);
        let rot_inv_vector = gyrovector_inv(rot_vector.clone(), k);

        let tmp_point = gyrovector_add(point.clone(), rot_vector, k);
        gyrovector_add(tmp_point, rot_inv_vector, k)
    }

    pub fn gyrovector_slerp(a: Point, b: Point, t: f64, k: f64) -> Point {
        let angle = gyrovector_distance(a, b, k);
        let sin_angle = angle.sinh();
    
        if sin_angle.abs() < f64::EPSILON {
            return gyrovector_lerp(a, b, t, k);
        }
    
        let coeff_a = ((1.0 - t) * angle).sinh() / sin_angle;
        let coeff_b = (t * angle).sinh() / sin_angle;
    
        let term_a = gyrovector_scalar_mul(coeff_a, a, k);
        let term_b = gyrovector_scalar_mul(coeff_b, b, k);
    
        gyrovector_add(term_a, term_b, k)
    }
}

mod H2xE {
    use crate::gyrovector;
    
    use nalgebra as na;

    pub struct H2xEPoint {
        pub coordinates: na::Vector3<f64>,
    }
    
    impl H2xEPoint {
        pub fn new(x: f64, y: f64, z: f64) -> Self {
            H2xEPoint {
                coordinates: na::Vector3::new(x, y, z),
            }
        }
    }
    
    pub fn combined_h2xe_distance(a: &H2xEPoint, b: &H2xEPoint, k: f64) -> f64 {
        let xy_distance = gyrovector::gyrovector_distance(
            na::Vector3::new(a.coordinates.x, a.coordinates.y, 0.0),
            na::Vector3::new(b.coordinates.x, b.coordinates.y, 0.0),
            k,
        );
        let z_distance = (b.coordinates.z - a.coordinates.z).abs();
    
        // Root of the sum of squares (RSS) of hyperbolic XY distance and Euclidean Z distance
        (xy_distance.powi(2) + z_distance.powi(2)).sqrt()
    }
    
    pub fn h2xe_interpolate(a: &H2xEPoint, b: &H2xEPoint, t: f64, k: f64) -> H2xEPoint {
        let xy_interp = gyrovector::gyrovector_slerp(
            na::Vector3::new(a.coordinates.x, a.coordinates.y, 0.0),
            na::Vector3::new(b.coordinates.x, b.coordinates.y, 0.0),
            t,
            k,
        );
        let z_interp = a.coordinates.z + t * (b.coordinates.z - a.coordinates.z);
    
        H2xEPoint::new(xy_interp.x, xy_interp.y, z_interp)
    }
}


#[cfg(test)]
mod tests {
    use crate::gyrovector;
    use crate::H2xE;

    use super::*;

    use nalgebra as na;

    #[test]
    fn it_works() {
        let k = -1.0;
        let point_a = H2xE::H2xEPoint::new(1.0, 2.0, 3.0);
        let point_b = H2xE::H2xEPoint::new(4.0, 6.0, 9.0);
        let expected_distance = 15.416214574631766; //TODO: actually write the test

        let calculated_distance = H2xE::combined_h2xe_distance(&point_a, &point_b, k);

        assert!(
            (calculated_distance - expected_distance).abs() < 1e-6,
            "Expected distance: {}, Calculated distance: {}",
            expected_distance,
            calculated_distance
        );
    }
}

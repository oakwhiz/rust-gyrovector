//This code has not yet been tested. It may not be correct.

mod gyrovector {
    use nalgebra as na;

    pub(crate) type Point = na::Vector3<f64>;

    pub fn gyrovector_add(u: Point, v: Point, k: f64) -> Point {
        let uv = u.dot(&v);
    
        if k.abs() < f64::EPSILON {
            return u + v;
        }
    
        let denominator = 1.0 - k * (2.0 * uv + uv * uv);
        u + (1.0 / k + uv) * ((1.0 / k + uv + uv * uv) / denominator) * (v - u)
    }
    
    // Gyrovector scalar multiplication with curvature constant K
    /*pub fn gyrovector_scalar_mul(scalar: f64, v: Point, k: f64) -> Point {
        if k.abs() < f64::EPSILON {
            return scalar * v;
        }
    
        let norm_v = v.norm();
        let lambda = (scalar * norm_v).tanh();
        (v / norm_v) * lambda
    }*/
    pub fn gyrovector_scalar_mul(scalar: f64, v: Point, k: f64) -> Point {
        if k.abs() < f64::EPSILON {
            return scalar * v;
        }
    
        let norm_v = v.norm();
        let lambda = (k * scalar * norm_v).tanh();
        (v / norm_v) * lambda
    }
    
    
    // Gyrovector scalar division with curvature constant K
    /*pub fn gyrovector_scalar_div(v: Point, scalar: f64, k: f64) -> Point {
        if k.abs() < f64::EPSILON {
            return v / scalar;
        }
    
        let norm_v = v.norm();
        let lambda = (norm_v / scalar).tanh();
        (v / norm_v) * lambda
    }*/
    pub fn gyrovector_scalar_div(v: Point, scalar: f64, k: f64) -> Point {
        if k.abs() < f64::EPSILON {
            return v / scalar;
        }
    
        let norm_v = v.norm();
        let lambda = (k * norm_v).tanh() / scalar;
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

    pub fn gyrovector_log(a: Point, b: Point, k: f64) -> Point {
        let inv_a = gyrovector_inv(a, k);
        let inv_a_b = gyrovector_add(inv_a, b, k);
        gyrovector_scalar_mul(1.0 / k, inv_a_b, k)
    }
    
    /*
    pub fn gyrovector_exp(a: Point, v: Point, k: f64) -> Point {
        let k_v = gyrovector_scalar_mul(k, v, k);
        let exp_k_v = gyrovector_add(a, k_v, k);
        gyrovector_scalar_mul(1.0 / k, exp_k_v, k)
    }
    */
    pub fn gyrovector_exp(a: Point, v: Point, k: f64) -> Point {
        let norm_v = v.norm();
        if norm_v < std::f64::EPSILON {
            return a;
        }
        let cosh_term = (k * norm_v).cosh();
        let sinh_term = (k * norm_v).sinh() / (k * norm_v);
        let new_v = v * sinh_term;
        let scale = (1.0 + k * a.norm_squared()) * cosh_term + a.dot(&v) * sinh_term;
        let result = (a + new_v) / scale;
        result
    }
    

    /*
    pub fn gyrovector_mul(a: Point, b: Point, k: f64) -> Point {
        let dot_product = a.dot(&b);
        let alpha = (k * dot_product).tanh();
        let sum = a + b;
        let scale = 1.0 / (1.0 + k * dot_product);
        sum * scale
    }
    */

    pub fn gyrovector_mul(a: Point, b: Point, k: f64) -> Point {
        let num = a.norm_squared() * b.norm_squared() - 2.0 * k * a.dot(&b) - k * k * (a.cross(&b)).norm_squared();
        let denom = 1.0 - 2.0 * k * a.dot(&b) + k * k * a.norm_squared() * b.norm_squared();
        let scale = num / denom;
        let cross_term = a.cross(&b);
        let sum = a + b + k * cross_term;
    
        Point::new(
            sum.x * scale,
            sum.y * scale,
            sum.z * scale,
        )
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
    /*
    pub fn gyrovector_lerp(a: Point, b: Point, t: f64, k: f64) -> Point {
        let d = gyrovector_distance(a, b, k);
        let s = t * d;
        let direction = gyrovector_scalar_mul(1.0 / d, gyrovector_add(b, gyrovector_inv(a, k), k), k);
        gyrovector_add(a, gyrovector_scalar_mul(s, direction, k), k)
    }
    */

    /*
    pub fn gyrovector_lerp(a: Point, b: Point, t: f64, k: f64) -> Point {
        let d = gyrovector_distance(a, b, k);
        let s = t * d;
        let direction = gyrovector_scalar_mul(1.0 / d, gyrovector_add(b, gyrovector_inv(a, k), k), k);
        let result = gyrovector_add(a, gyrovector_scalar_mul(s, direction, k), k);
    
        println!("a: {:?}", a);
        println!("b: {:?}", b);
        println!("d: {}", d);
        println!("s: {}", s);
        println!("direction: {:?}", direction);
        println!("result: {:?}", result);
    
        result
    }
    */

    /*
    pub fn gyrovector_lerp(a: Point, b: Point, t: f64, k: f64) -> Point {
        let d = gyrovector_distance(a, b, k);
        let s = t * d;
        let direction = gyrovector_scalar_mul(1.0 / d, gyrovector_add(b, gyrovector_inv(a, k), k), k);
        gyrovector_add(a, gyrovector_scalar_mul(s / k, direction, k), k)
    }*/

    /*
    pub fn gyrovector_lerp(a: Point, b: Point, t: f64, k: f64) -> Point {
        let log_ab = gyrovector_log(a, b, k);
        let t_log_ab = gyrovector_scalar_mul(t, log_ab, k);
        gyrovector_exp(a, t_log_ab, k)
    }
    */
    
    pub fn gyrovector_lerp(a: Point, b: Point, t: f64, k: f64) -> Point {
        let log_ba = gyrovector_log(a, b, k);
        let t_log_ba = gyrovector_scalar_mul(t, log_ba, k);
        gyrovector_mul(a, gyrovector_exp(a, t_log_ba, k), k)
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

    use assert_approx_eq::assert_approx_eq;

    #[test]
    fn it_works() {
        //expected results on this test has not yet been checked
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

    /*
    #[test]
    fn test_gyrovector_lerp_midpoint() {
        let k = -1.0;

        let a = gyrovector::Point::new(-1.0, 1.0, 0.0);
        let b = gyrovector::Point::new(1.0, 1.0, 0.0);
        let c = gyrovector::Point::new(0.0, -1.0, 0.0);

        let m_expected = gyrovector::Point::new(0.5, 0.0, 0.0);

        let m_computed = gyrovector::gyrovector_lerp(b, c, 0.5, k);

        assert_approx_eq!(m_expected.x, m_computed.x, 1e-6);
        assert_approx_eq!(m_expected.y, m_computed.y, 1e-6);
        assert_approx_eq!(m_expected.z, m_computed.z, 1e-6);
    }
    */

    /*
    #[test]
    fn test_gyrovector_lerp_midpoint() {
        let a = gyrovector::Point::new(1.0, 135.0_f64.to_radians());
        let b = gyrovector::Point::new(1.0, 45.0_f64.to_radians());
        let c = gyrovector::Point::new(1.0, 270.0_f64.to_radians());

        let k = -1.0;

        let midpoint_bc = gyrovector::gyrovector_lerp(b, c, 0.5, k);

        // Check if the midpoint of BC is equidistant from B and C
        let distance_mb = H2xE::combined_h2xe_distance(midpoint_bc, b);
        let distance_mc = H2xE::combined_h2xe_distance(midpoint_bc, c);

        let diff = (distance_mb - distance_mc).abs();
        assert!((diff - 1e-6).abs() < 1e-6, "assertion failed: `(left !== right)` (left: `{:?}`, right: `{:?}`, expect diff: `1e-6`, real diff: `{:?}`)", distance_mb, distance_mc, diff);
    }
    */

    /*
    #[test]
    fn test_gyrovector_lerp_midpoint() {
        let a = gyrovector::Point::new(1.0, 135.0_f64.to_radians(), 0.0);
        let b = gyrovector::Point::new(1.0, 45.0_f64.to_radians(), 0.0);
        let c = gyrovector::Point::new(1.0, 270.0_f64.to_radians(), 0.0);

        let k = -1.0;

        let midpoint_bc = gyrovector::gyrovector_lerp(b, c, 0.5, k);

        // Check if the midpoint of BC is equidistant from B and C
        let distance_mb = gyrovector::gyrovector_distance(midpoint_bc, b, k);
        let distance_mc = gyrovector::gyrovector_distance(midpoint_bc, c, k);

        let diff = (distance_mb - distance_mc).abs();
        assert!((diff - 1e-6).abs() < 1e-6, "assertion failed: `(left !== right)` (left: `{:?}`, right: `{:?}`, expect diff: `1e-6`, real diff: `{:?}`)", distance_mb, distance_mc, diff);
    }
    */

    /*
    #[test]
    fn test_gyrovector_lerp_midpoint() {
        let a = gyrovector::Point::new(1.0, 1.0, 0.0);
        let b = gyrovector::Point::new(3.0, 1.0, 0.0);
        let c = gyrovector::Point::new(3.0, 3.0, 0.0);

        let k = -1.0;

        let midpoint_bc = gyrovector::gyrovector_lerp(b, c, 0.5, k);

        // Expected midpoint of BC
        let m_expected = gyrovector::Point::new(3.0, 1.927050983124842, 0.0);

        assert_approx_eq!(m_expected.x, midpoint_bc.x, 1e-6);
        assert_approx_eq!(m_expected.y, midpoint_bc.y, 1e-6);
        assert_approx_eq!(m_expected.z, midpoint_bc.z, 1e-6);
    }
    */

    /*
    #[test]
    fn test_gyrovector_lerp_midpoint() {
        let a = gyrovector::Point::new(1.0, 135.0_f64.to_radians(), 0.0);
        let b = gyrovector::Point::new(1.0, 45.0_f64.to_radians(), 0.0);
        let c = gyrovector::Point::new(1.0, 270.0_f64.to_radians(), 0.0);

        let k = -1.0;

        let midpoint_bc = gyrovector::gyrovector_lerp(b, c, 0.5, k);

        // Update the expected values
        let m_expected = gyrovector::Point::new(3.00019, 2.00012, 0.0);

        assert_approx_eq!(m_expected.x, midpoint_bc.x, 1e-6);
        assert_approx_eq!(m_expected.y, midpoint_bc.y, 1e-6);
        assert_approx_eq!(m_expected.z, midpoint_bc.z, 1e-6);
    }
    */

    /*

    #[test]
    fn test_gyrovector_lerp_midpoint() {
        let a = gyrovector::Point::new(1.0, 1.0, 0.0);
        let b = gyrovector::Point::new(5.0, 3.0, 0.0);
        let t = 0.5;
        let k = -1.0;

        let result = gyrovector::gyrovector_lerp(a, b, t, k);

        println!("a: {:?}", a);
        println!("b: {:?}", b);
        println!("t: {}", t);
        println!("k: {}", k);
        println!("result: {:?}", result);

        assert_approx_eq!(result.x, 3.00019, 1e-6);
        assert_approx_eq!(result.y, 2.00012, 1e-6);
    }
    */

    #[test]
    fn test_gyrovector_lerp_midpoint() {
        let a = gyrovector::Point::new(1.0, 1.0, 0.0);
        let b = gyrovector::Point::new(5.0, 3.0, 0.0);
        let t = 0.5;
        let k = -1.0;

        let result = gyrovector::gyrovector_lerp(a, b, t, k);

        println!("a: {:?}", a);
        println!("b: {:?}", b);
        println!("t: {}", t);
        println!("k: {}", k);
        println!("result: {:?}", result);

        assert_approx_eq!(result.x, 2.9999999999999996, 1e-6);
        assert_approx_eq!(result.y, 2.0, 1e-6);
    }




}

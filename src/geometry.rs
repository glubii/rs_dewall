use nalgebra::distance_squared;
use nalgebra::Point2;
use std::fmt::Debug;
use std::ops::{Deref, DerefMut};

#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Point(Point2<f32>);

impl Point {
    pub fn new(x: f32, y: f32) -> Self {
        Point(Point2::new(x, y))
    }

    pub fn distance_squared(&self, other: &Point) -> f32 {
        distance_squared(&self.0, &other.0)
    }

    //https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Line_defined_by_two_points
    //SIGNED, use abs for min !
    pub fn distance_to_edge(&self, edge: &Edge) -> f32 {
        let f = (edge.end.x - edge.start.x) * (edge.start.y - self.y)
            - (edge.start.x - self.x) * (edge.end.y - edge.start.y);

        let g =
            ((edge.end.x - edge.start.x).powf(2.0) + (edge.end.y - edge.start.y).powf(2.0)).sqrt();

        return f / g;
    }

    /**
     * @brief Checks if the point is on the left side of the edge.
     * Done by checking if the cross prod. / det of the point matrix > 0
     *
     * @returns True if on the left side
     */
    pub fn left_side_of_edge(&self, edge: &Edge) -> bool {
        ((edge.end.x - edge.start.x) * (self.y - edge.start.y)
            - (edge.end.y - edge.start.y) * (self.x - edge.start.x))
            >= 0.0
    }

    /*
     * @brief Returns the center and the radius of the circumcircle around the three points
     *
     * URL (version: 2022-02-22): https://math.stackexchange.com/q/3503338
     * Scott (https://math.stackexchange.com/users/740203/scott)
     * https://www.geeksforgeeks.org/equation-of-circle-when-three-points-on-the-circle-are-given/
     */
    pub fn circumcircle(&self, p2: &Point, p3: &Point) -> (Point, f32) {
        if self == p2 || p2 == p3 || p3 == self {
            panic!("The points may not be the same")
        }

        let x12 = self.x - p2.x;
        let x13 = self.x - p3.x;

        let y12 = self.y - p2.y;
        let y13 = self.y - p3.y;

        let x31 = p3.x - self.x;
        let x21 = p2.x - self.x;

        let y31 = p3.y - self.y;
        let y21 = p2.y - self.y;

        // x1^2 - x3^2
        let sx13 = self.x.powf(2.0) - p3.x.powf(2.0);
        let sy13 = self.y.powf(2.0) - p3.y.powf(2.0);
        let sx21 = p2.x.powf(2.0) - self.x.powf(2.0);
        let sy21 = p2.y.powf(2.0) - self.y.powf(2.0);

        let f = ((sx13) * (x12) + (sy13) * (x12) + (sx21) * (x13) + (sy21) * (x13))
            / (2.0 * ((y31) * (x12) - (y21) * (x13)));
        let g = ((sx13) * (y12) + (sy13) * (y12) + (sx21) * (y13) + (sy21) * (y13))
            / (2.0 * ((x31) * (y12) - (x21) * (y13)));

        let c = -(self.x.powf(2.0)) - (self.y.powf(2.0)) - 2.0 * g * (self.x) - 2.0 * f * (self.y);

        // eqn of circle be x^2 + y^2 + 2*g*x + 2*f*y + c = 0
        // where centre is (h = -g, k = -f) and radius r
        // as r^2 = h^2 + k^2 - c
        let h = -g;
        let k = -f;
        let sqr_of_r = h * h + k * k - c;

        // r is the radius
        let r = sqr_of_r.sqrt();

        return (Point::new(h, k), r);
    }
}

impl Deref for Point {
    type Target = Point2<f32>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Point {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Edge {
    pub start: Point,
    pub end: Point,
}

impl Edge {
    pub fn new(start: Point, end: Point) -> Self {
        if start == end {
            panic!("Start- and endpoint of an edge may not be the same");
        }
        Edge { start, end }
    }
}

//Make the comparison direction invariant
impl PartialEq for Edge {
    fn eq(&self, other: &Self) -> bool {
        self.start == other.start && self.end == other.end
            || self.end == other.start && self.start == other.end
    }
}

#[derive(Clone, Copy)]
pub struct Triangle {
    pub p1: Point,
    pub p2: Point,
    pub p3: Point,

    pub e1: Edge,
    pub e2: Edge,
    pub e3: Edge,
}

impl Triangle {
    pub fn new(p1: Point, p2: Point, p3: Point) -> Self {
        if p1 == p2 || p2 == p3 || p3 == p1 {
            panic!("Edges of a triangle may not be the same points");
        }
        let e1 = Edge::new(p2, p3);
        let e2 = Edge::new(p3, p1);
        let e3 = Edge::new(p1, p2);

        assert!(p1.left_side_of_edge(&e1));
        assert!(p2.left_side_of_edge(&e2));
        assert!(p3.left_side_of_edge(&e3));

        Triangle {
            p1,
            p2,
            p3,
            e1,
            e2,
            e3,
        }
    }
    pub fn point_inside(&self, point: Point) -> bool {
        let e1 = Edge::new(self.p1, self.p2);
        let e2 = Edge::new(self.p2, self.p3);
        let e3 = Edge::new(self.p3, self.p1);
        return point.left_side_of_edge(&e1)
            && point.left_side_of_edge(&e2)
            && point.left_side_of_edge(&e3);
    }
}

pub struct SpacePartition {
    pub span_x: Edge,
    pub span_y: Edge,
    pub alpha: Edge,
}

impl SpacePartition {
    pub fn new(span_x: Edge, span_y: Edge, alpha: Edge) -> Self {
        SpacePartition {
            span_x,
            span_y,
            alpha,
        }
    }

    pub fn partition(&self) -> (SpacePartition, SpacePartition) {
        let alpha1: Edge;
        let alpha2: Edge;
        let span_x_1: Edge;
        let span_x_2: Edge;
        let span_y_1: Edge;
        let span_y_2: Edge;

        if self.alpha.start.x == self.alpha.end.x {
            let new_y = self.alpha.start.y + (self.alpha.end.y - self.alpha.start.y) / 2.0;
            alpha2 = Edge::new(
                Point::new(self.span_x.start.x, new_y),
                Point::new(
                    self.span_x.start.x + (self.span_x.end.x - self.span_x.start.x) / 2.0,
                    new_y,
                ),
            );
            alpha1 = Edge::new(alpha2.end, Point::new(self.span_x.end.x, new_y));
            span_x_2 = Edge::new(
                self.span_x.start,
                Point::new(
                    self.span_x.start.x + (self.span_x.end.x - self.span_x.start.x) / 2.0,
                    self.span_x.start.y,
                ),
            );
            span_x_1 = Edge::new(span_x_2.end, self.span_x.end);

            span_y_1 = self.span_y;
            span_y_2 = self.span_y;
        } else {
            let new_x = self.alpha.start.x + (self.alpha.end.x - self.alpha.start.x) / 2.0;
            alpha1 = Edge::new(
                Point::new(new_x, self.span_y.start.y),
                Point::new(
                    new_x,
                    self.span_y.start.y + (self.span_y.end.y - self.span_y.start.y) / 2.0,
                ),
            );
            alpha2 = Edge::new(alpha1.end, Point::new(new_x, self.span_y.end.y));

            span_y_1 = Edge::new(
                self.span_y.start,
                Point::new(
                    self.span_y.start.x,
                    self.span_y.start.y + (self.span_y.end.y - self.span_y.start.y) / 2.0,
                ),
            );
            span_y_2 = Edge::new(span_y_1.end, self.span_y.end);

            span_x_1 = self.span_x;
            span_x_2 = self.span_x;
        }

        return (
            SpacePartition::new(span_x_1, span_y_1, alpha1),
            SpacePartition::new(span_x_2, span_y_2, alpha2),
        );
    }
}

#[cfg(test)]
mod tests {
    use crate::Edge;
    use crate::Point;
    use crate::SpacePartition;

    #[test]
    fn test_length() {
        let p1 = Point::new(1.0, 1.0);
        let p2 = Point::new(3.0, 1.0);

        assert_eq!(p1.distance(&p2), 2.0);
    }

    #[test]
    fn test_distance_to_edge() {
        let p1 = Point::new(1.0, 1.0);
        let p2 = Point::new(3.0, 1.0);
        let p3 = Point::new(2.0, 2.0);

        let edge = Edge::new(p1, p2);
        let inv_edge = Edge::new(p2, p1);

        assert_eq!(p3.distance_to_edge(&edge), -1.0);
        assert_eq!(p3.distance_to_edge(&inv_edge), 1.0);
    }

    #[test]
    fn test_left_side_of_edge() {
        let p1 = Point::new(1.0, 1.0);
        let p2 = Point::new(3.0, 1.0);
        let p3 = Point::new(2.0, 2.0);

        let edge = Edge::new(p1, p2);
        let inv_edge = Edge::new(p2, p1);

        assert_eq!(p3.left_side_of_edge(&edge), true);
        assert_eq!(p3.left_side_of_edge(&inv_edge), false);
    }

    #[test]
    fn test_circumcircle() {
        let p1 = Point::new(1.0, 1.0);
        let p2 = Point::new(3.0, 1.0);
        let p3 = Point::new(2.0, 2.0);
        let center = Point::new(2.0, 1.0);
        let (c, r) = p1.circumcircle(&p2, &p3);

        assert_eq!(r, 1.0);
        assert_eq!(c, center);
    }

    #[test]
    fn test_space_partition_x() {
        let p1 = Point::new(300.0, 100.0);
        let p2 = Point::new(700.0, 100.0);
        let p3 = Point::new(300.0, 300.0);

        let span_x = Edge::new(p1, p2);
        let span_y = Edge::new(p1, p3);

        let alpha = Edge::new(Point::new(500.0, 100.0), Point::new(500.0, 300.0));

        let alpha_2 = Edge::new(Point::new(300.0, 200.0), Point::new(500.0, 200.0));
        let alpha_1 = Edge::new(Point::new(500.0, 200.0), Point::new(700.0, 200.0));

        let span_x_2 = Edge::new(p1, alpha.start);
        let span_x_1 = Edge::new(alpha.start, p2);

        let sp = SpacePartition::new(span_x, span_y, alpha);
        let (sp1, sp2) = sp.partition();
        assert_eq!(sp1.alpha, alpha_1);
        assert_eq!(sp2.alpha, alpha_2);

        assert_eq!(sp1.span_x, span_x_1);
        assert_eq!(sp2.span_x, span_x_2);

        assert_eq!(sp1.span_y, span_y);
        assert_eq!(sp2.span_y, span_y);
    }

    #[test]
    fn test_space_partition_y() {
        let p1 = Point::new(100.0, 300.0);
        let p2 = Point::new(100.0, 700.0);
        let p3 = Point::new(300.0, 300.0);

        let span_y = Edge::new(p1, p2);
        let span_x = Edge::new(p1, p3);

        let alpha = Edge::new(Point::new(100.0, 500.0), Point::new(300.0, 500.0));
        let alpha_1 = Edge::new(Point::new(200.0, 300.0), Point::new(200.0, 500.0));
        let alpha_2 = Edge::new(Point::new(200.0, 500.0), Point::new(200.0, 700.0));

        let span_y_1 = Edge::new(p1, alpha.start);
        let span_y_2 = Edge::new(alpha.start, p2);

        let sp = SpacePartition::new(span_x, span_y, alpha);
        let (sp1, sp2) = sp.partition();
        assert_eq!(sp1.alpha, alpha_1);
        assert_eq!(sp2.alpha, alpha_2);

        assert_eq!(sp1.span_x, span_x);
        assert_eq!(sp2.span_x, span_x);

        assert_eq!(sp1.span_y, span_y_1);
        assert_eq!(sp2.span_y, span_y_2);
    }
}

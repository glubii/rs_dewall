
#[derive(Clone, Copy, PartialEq, Eq, Ord)]
pub struct Point {
    pub x: u32,
    pub y: u32
}

impl Point {
    pub fn new (x: u32, y: u32) -> Self {
        Point{x, y}
    }

    pub fn length(&self, p2: &Point) -> f32 {
        ((p2.x as f32 - self.x as f32).powf(2.0) + (p2.y as f32 - self.y as f32).powf(2.0)).sqrt()
    }

    pub fn distance_no_sqrt(&self, p2: &Point) -> f32 {
        (p2.x as f32 - self.x as f32).powf(2.0) + (p2.y as f32 - self.y as f32).powf(2.0)
    }

    //https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Line_defined_by_two_points
    //SIGNED, use abs for min !
    pub fn distance_to_edge(&self, edge: &Edge) -> f32 {
        let f = (
            (edge.end.x as f32 - edge.start.x as f32) *(edge.start.y as f32 - self.y as f32)
          - (edge.start.x as f32 - self.x as f32) * (edge.end.y as f32 - edge.start.y as f32)
        );

        let g = (
            (edge.end.x as f32 - edge.start.x as f32).powf(2.0)
          + (edge.end.y as f32 - edge.start.y as f32).powf(2.0)
        ).sqrt();

        return f/g;
    }
    
    //Checks if the cross prod. / det of the point matrix > 0
    pub fn left_side_of_edge(&self, edge: &Edge) -> bool {
        ( (edge.end.x as f32 - edge.start.x as f32)
        *(self.y  as f32- edge.start.y as f32)
        - (edge.end.y as f32 - edge.start.y as f32)
        *(self.x as f32 - edge.start.x as f32) ) > 0.0
    }

    /*
     * Get the equation of a circle when given 3 points
     * URL (version: 2022-02-22): https://math.stackexchange.com/q/3503338
     * Scott (https://math.stackexchange.com/users/740203/scott)
     * 
     * Returns the center and the radius of the circumcircle
     */

    pub fn circumcircle(&self, p2: &Point, p3: &Point) -> (Point, f32) {
        if self == p2 || p2 == p3 || p3 == self {
            //Error, return max radius
            return (Point::new(0, 0), 999.9)
        }

        //https://www.geeksforgeeks.org/equation-of-circle-when-three-points-on-the-circle-are-given/

        let x12 = self.x as f32 - p2.x as f32;
        let x13 = self.x as f32 - p3.x as f32;

        let y12 = self.y as f32 - p2.y as f32;
        let y13 = self.y as f32 - p3.y as f32;
        
        let x31 = p3.x as f32 - self.x as f32;
        let x21 = p2.x as f32 - self.x as f32;

        let y31 = p3.y as f32 - self.y as f32;
        let y21 = p2.y as f32 - self.y as f32;

        // x1^2 - x3^2
        let sx13 = self.x.pow(2) as f32 - p3.x.pow(2) as f32;
        let sy13 = self.y.pow(2) as f32 - p3.y.pow(2) as f32;
        let sx21 = p2.x.pow(2) as f32 - self.x.pow(2) as f32;
        let sy21 = p2.y.pow(2) as f32 - self.y.pow(2) as f32;
 
        let f = ((sx13) * (x12)
             + (sy13) * (x12)
             + (sx21) * (x13)
             + (sy21) * (x13))
            / (2.0 * ((y31) * (x12) - (y21) * (x13)));
        let g = ((sx13) * (y12)
             + (sy13) * (y12)
             + (sx21) * (y13)
             + (sy21) * (y13))
            / (2.0 * ((x31) * (y12) - (x21) * (y13)));
 
        let c = -(self.x.pow(2) as f32) - (self.y.pow(2) as f32) - 2.0 * g * (self.x as f32) - 2.0 * f * (self.y as f32);
 
        // eqn of circle be x^2 + y^2 + 2*g*x + 2*f*y + c = 0
        // where centre is (h = -g, k = -f) and radius r
        // as r^2 = h^2 + k^2 - c
        let h = -g;
        let k = -f;
        let sqr_of_r = h * h + k * k - c;
 
        // r is the radius
        let r = sqr_of_r.sqrt();
        
        return (Point::new(h as u32, k as u32), r)
    }
}

impl PartialOrd for Point {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
       if self.x == other.x {
           return Some(self.y.cmp(&other.y))
       }
       Some(self.x.cmp(&other.x))
    }
}

#[derive(Clone, Copy)]
pub struct Edge {
    pub start: Point,
    pub end: Point
}

impl Edge {
    pub fn new (start: Point, end: Point) -> Self {
        Edge{start, end}
    }

    //TODO
    //fn len()
}

//Make the comparison direction invariant
impl PartialEq for Edge {
    fn eq(&self, other: &Self) -> bool {
        self.start == other.start && self.end == other.end ||
        self.end == other.start && self.start == other.end 
    }
}

#[derive(Clone, Copy)]
pub struct Triangle{
    pub p1: Point,
    pub p2: Point,
    pub p3: Point
}

impl Triangle {
    pub fn new (p1: Point, p2: Point, p3: Point) -> Self {
        Triangle{p1, p2, p3}
    }
    pub fn point_inside(&self, point: Point) -> bool {
        let e1 = Edge::new(self.p1, self.p2);
        let e2 = Edge::new(self.p2, self.p3);
        let e3 = Edge::new(self.p3, self.p1);
        return point.left_side_of_edge(&e1) && point.left_side_of_edge(&e2) && point.left_side_of_edge(&e3);
    }
}


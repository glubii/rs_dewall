use std::mem::align_of_val;

use image;
use rand::Rng;

#[derive(Clone, Copy, PartialEq, Eq, Ord)]
struct Point {
    pub x: u32,
    pub y: u32
}

impl Point {
    fn new (x: u32, y: u32) -> Self {
        Point{x, y}
    }

    fn length(&self, p2: &Point) -> f32 {
        (p2.x as f32 - self.x as f32).powf(2.0) + (p2.y as f32 - self.y as f32).powf(2.0).sqrt()
    }

    fn distance_no_sqrt(&self, p2: &Point) -> f32 {
        (p2.x as f32 - self.x as f32).powf(2.0) + (p2.y as f32 - self.y as f32).powf(2.0)
    }

    //TODO
    //fn dist to edge

    fn radius_min_circle(&self, p2: &Point, p3: &Point) -> f32 {
        //https://math.stackexchange.com/questions/213658/get-the-equation-of-a-circle-when-given-3-points
        p2.length(p3) * self.length(p2) * self.length(p3)
        / (2.0 * (self.x as i32 * p2.y as i32 + p2.x as i32 * p3.y as i32 + p3.x as i32 * self.y as i32 - self.y as i32 * p2.x as i32 - p2.y as i32 * p3.x as i32 - p3.y as i32 * self.x as i32).abs() as f32)
        //TODO -> Triangle.area
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
struct Edge {
    pub start: Point,
    pub end: Point
}

impl Edge {
    fn new (start: Point, end: Point) -> Self {
        Edge{start, end}
    }

    //TODO
    //fn len()
}

#[derive(Clone, Copy)]
struct Triangle{
    pub p1: Point,
    pub p2: Point,
    pub p3: Point
}

impl Triangle {
    fn new (p1: Point, p2: Point, p3: Point) -> Self {
        Triangle{p1, p2, p3}
    }
}

//See Cignoni et al
fn de_wall(points: Vec<Point>, wall: Edge, afl: Vec<Triangle>) -> Vec<Triangle> {
    //https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Line_defined_by_two_points
    let distance = |p: &Point|
        ((wall.end.x as f32 - wall.start.x as f32) *(wall.start.y as f32 - p.y as f32)
        - (wall.start.x as f32 - p.x as f32) * (wall.end.y as f32 - wall.start.y as f32))
        / ((wall.end.x as f32- wall.start.x as f32).powf(2.0)
        + (wall.end.y as f32 - wall.start.y as f32).powf(2.0)).sqrt();
 
    //Divide Points into lower (negative) and upper (positive) groups
    let mut upper_points: Vec<Point> = Vec::new();
    let mut lower_points: Vec<Point> = Vec::new();
    for p in &points {
        if distance(p) > 0.0 {
            upper_points.push(*p);
        }
        else {
            lower_points.push(*p);
        }
    }

    let mut d_simplex_list: Vec<Triangle> = Vec::new();
    let mut afl: Vec<Edge> = Vec::new();
    let mut afl_sigma: Vec<Edge> = Vec::new();
    let mut afl_1: Vec<Edge> = Vec::new();
    let mut afl_2: Vec<Edge> = Vec::new();

    if afl.is_empty() {
        //Make first Simplex by finding the point closest to the intersecting plane and the nearest one on the opposite side
      //  let minimum = |p1: &Point, p2: &Point| ( if distance(p1) < distance(p2) {p1} else {p2} );

        //points.iter().reduce(minimum);
        let point_first_opt = points.iter().reduce(|p1, p2| ( if distance(p1).abs() < distance(p2).abs() {p1} else {p2}));

        let point_first = match point_first_opt {
            Some(x) => {x}
            None => {panic!()} //Idk if this should ever occur
        };

        println!("{} {}", point_first.x, point_first.y); 

        //Figure out via the cross product on which side of the intersection the point is
        let point_second_opt:Option<&Point>;
        if distance(point_first) > 0.0 {
            point_second_opt = lower_points.iter().reduce(
                |p1, p2|
                (
                    if p1.distance_no_sqrt(point_first) < p2.distance_no_sqrt(point_first)
                        {p1}
                    else
                        {p2}
                )
            );
        }
        else {
            point_second_opt = upper_points.iter().reduce(
                |p1, p2|
                (
                    if p1.distance_no_sqrt(point_first) < p2.distance_no_sqrt(point_first)
                        {p1}
                    else
                        {p2}
                )
            );
        }
        
        let point_second= match point_second_opt{
            Some(x) => {x}
            None => {panic!()} //Idk if this should ever occur
        };

        println!("{} {}", point_second.x, point_second.y); 

        let point_third_opt = points.iter().reduce(
            |p1, p2|
            (
                if p2.radius_min_circle(point_first, point_second) < p1.radius_min_circle(point_first, point_second) {
        
                    if p2 != point_first && p2 != point_second {
                        p2
                    }
                    else {
                        p1
                    }
                }
                else {
                    p1
                }
            ));
        
        let point_third = match point_third_opt{
            Some(x) => {x}
            None => {panic!()} //Idk if this should ever occur
        };
        
        println!("{} {}", point_third.x, point_third.y);

        afl.push(Edge::new(*point_first, *point_second));
        afl.push(Edge::new(*point_second, *point_third));
        afl.push(Edge::new(*point_third, *point_first));
        d_simplex_list.push(Triangle::new(*point_first, *point_second, *point_third));
    }

    for face in afl {
        if lower_points.contains(&face.start) && lower_points.contains(&face.end){
            afl_1.push(face);
        }

        if upper_points.contains(&face.start) && upper_points.contains(&face.end){
            afl_2.push(face);
        }
        else{
            //splitting wall intersects triangle
            afl_sigma.push(face);
        }
    }

    //for face in afl_sigma {
    while !afl_sigma.is_empty() {
        let face_opt = afl_sigma.pop();
        
        let face = match face_opt{
            Some(x) => {x}
            None => {panic!()} //Idk if this should ever occur
        };

        //let t_opt: Option<Triangle> = simplex(face, points);


        let point_third_opt = points.iter().reduce(
            |p1, p2|
            (
                if p2.radius_min_circle(&face.start, &face.end) < p1.radius_min_circle(&face.start, &face.end) {
        
                    if p2 != &face.start && p2 != &face.end{
                        p2
                    }
                    else {
                        p1
                    }
                }
                else {
                    p1
                }
            ));
        
        let point_third = match point_third_opt{
            Some(x) => {x}
            None => {panic!()} //Idk if this should ever occur
        };
        
        let triangle = Triangle::new(face.start, face.end, *point_third); 

        println!("{} {}", face.start.x, face.start.y);
        println!("{} {}", face.end.x, face.end.y);
        println!("{} {}", point_third.x, point_third.y);

        d_simplex_list.push(triangle);
        for new_face in [
            Edge::new(triangle.p1, triangle.p2),
            Edge::new(triangle.p2, triangle.p3),
            Edge::new(triangle.p3, triangle.p1)
        ] {
            if lower_points.contains(&new_face.start) && lower_points.contains(&new_face.end){
                afl_1.push(new_face);
            }
            if upper_points.contains(&new_face.start) && upper_points.contains(&new_face.end){
                afl_2.push(new_face);
            }
            else{
                //splitting wall intersects triangle
                afl_sigma.push(new_face);
            }
        }
    }



    //TODO: New walls
    //d_simplex_list.append(&mut de_wall(upper_points, wall, afl_1));
    //d_simplex_list.append(&mut de_wall(lower_points, wall, afl_2));

    return d_simplex_list;

}

fn main() {
    println!("Hello, world!");

    let imgx = 800;
    let imgy = 600;

    let mut imgbuf = image::RgbImage::new(imgx, imgy);

    for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
        let r = (255 * x/imgx) as u8;
        let g = (255 * y/imgy) as u8;
        
        *pixel = image::Rgb([r,g,0]);
    }

 
    //Random points
    let amount = 500;
    let mut rng = rand::thread_rng();

    let mut points_list: Vec<Point> = Vec::new();

    for _i in 0..amount {
        //TODO: Exclude the border region
        let x = rng.gen_range(0..imgx);
        let y = rng.gen_range(0..imgy);

        points_list.push(Point::new(x, y));

        let pixel = imgbuf.get_pixel_mut(x, y);
        *pixel = image::Rgb([255,255,255]);
    }

    //Corner Points
    
    //Sort
    points_list.sort();

    //Delauny
    let d_simplex_list = de_wall(points_list, Edge { start: Point::new(0, imgy/2), end: Point::new(imgx, imgy/2) },Vec::new());

    for tiangle in d_simplex_list {
        for point in [tiangle.p1, tiangle.p2, tiangle.p3] {
            let pixel = imgbuf.get_pixel_mut(point.x, point.y);
            *pixel = image::Rgb([0,0,255]);
        }
    }

    //colourize

    //export

    imgbuf.save("out.png").unwrap();
}

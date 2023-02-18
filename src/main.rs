use rand::Rng;

mod geometry;
use geometry::Point as Point;
use geometry::Edge as Edge;
use geometry::Triangle as Triangle;

/*
 * Fast divide and conquer Delaunay triangulation algorithm in E^d with d=2
 * Cignoni et al. 1998
 * 
 * P        the list of points to triangulate
 * x        TODO
 * y        TODO
 * alpha    A line which splits P in two disjoint Pointsets P1, P2
 * afl      TODO
 */
fn de_wall(P: Vec<Point>, x:Edge, y:Edge, alpha: Edge, afl: &mut Vec<Edge>) -> Vec<Triangle> {
    // (d-1) faces intersected by plane alpha
    let mut afl_sigma: Vec<Edge> = Vec::new();
    
    // (d-1) faces with all of the vertices in P1
    let mut afl_1: Vec<Edge> = Vec::new();
    
    // (d-1) faces with all of the vertices in P2
    let mut afl_2: Vec<Edge> = Vec::new();

    let mut P1: Vec<Point> = Vec::new();
    let mut P2: Vec<Point> = Vec::new();

    let mut sigma: Vec<Triangle> = Vec::new();

    //Divide Points into lower (negative) and upper (positive) groups
    //TODO: Use references instead?
    pointset_partition(&P, &alpha, &mut P1, &mut P2);

    println!("len of P {} {} {}", P.len(), P1.len(), P2.len());

    //Simplex Wall Construction
    if afl.is_empty() {
        let t = make_first_simplex(&P, &alpha, &P1, &P2);

        afl.push(Edge::new(t.p1, t.p2));
        afl.push(Edge::new(t.p2, t.p3));
        afl.push(Edge::new(t.p3, t.p1));
        sigma.push(t);
   
        
        println!("First simplex done:");
        println!("1 {} {}", t.p1.x, t.p1.y);
        println!("2 {} {}", t.p2.x, t.p2.y);
        println!("3 {} {}", t.p3.x, t.p3.y);
    }

    for f in afl {
        if P1.contains(&f.start) && P1.contains(&f.end){
            afl_1.push(*f);
        }

        else if P2.contains(&f.start) && P2.contains(&f.end){
            afl_2.push(*f);
        }
        else{
            //splitting wall intersects triangle
            afl_sigma.push(*f);
        }
    }

    //for _face in 0..6 {
    while !afl_sigma.is_empty() {

        // f = Extract(afl_sigma)
        println!("Loop start {}", afl_sigma.len());
        let face_opt = afl_sigma.pop();
        let f = match face_opt{
            Some(x) => {x}
            None => {panic!()} //Idk if this can ever occur
        };

        let t = make_simplex(&f, &P);

        /*
        println!("1 {} {}", t.p1.x, t.p1.y);
        println!("2 {} {}", t.p2.x, t.p2.y);
        println!("3 {} {}", t.p3.x, t.p3.y);
        */

        if true{ //t != null { TODO: Add option or something to make_simplex
            sigma.push(t);

            for f_apostroph in [
                Edge::new(t.p1,t.p2),
                Edge::new(t.p2,t.p3),
                Edge::new(t.p3,t.p1),
                ] {

                if f_apostroph != f {
                    if P1.contains(&f_apostroph.start) && P1.contains(&f_apostroph.end){
                        //println!("AFL1");
                        //println!("{}", afl_1.len());
                        update(&f_apostroph, &mut afl_1);
                        //println!("{}", afl_1.len());
                    }
                    else if P2.contains(&f_apostroph.start) && P2.contains(&f_apostroph.end){
                        //println!("AFL2");
                        //println!("{}", afl_2.len());
                        update(&f_apostroph, &mut afl_2);
                        //println!("{}", afl_2.len());
                    }
                    else{
                        //splitting wall intersects triangle
                        //println!("Sigma");
                        //println!("{}", afl_sigma.len());
                        update(&f_apostroph, &mut afl_sigma);
                        //println!("{}", afl_sigma.len());
                    }
                }
                else {
                    //println!("Eq")
                }
            }
        }
    }

    //TODO: Calculate new splitting planes
    //alpha1 = 

    //sigma.append(&mut de_wall(P1, alpha_1, afl_1));
    //sigma.append(&mut de_wall(P2, alpha_2, afl_2));

    return sigma;

}

//p.336
fn update (f: &Edge, L: &mut Vec<Edge>) {
    if L.contains(f) {
        L.retain(|&x| x != *f);
    }
    else {
        L.push(*f);
    }
}

//p.335
fn make_first_simplex(P: &Vec<Point>, alpha: &Edge, P1: &Vec<Point>, P2: &Vec<Point>) -> Triangle {
    //Find the point closest to alpha
    //let min_distance_to_alpha = |p1: &Point, p2: &Point| ( if p1.distance_to_edge(alpha) < p2.distance_to_edge(alpha) {p1} else {p2} );

    let corner_1 = match P.iter().reduce(
        |p1: &Point, p2: &Point| ( if p1.distance_to_edge(alpha).abs() < p2.distance_to_edge(alpha).abs() {p1} else {p2} )
    ) {
    //let p_1= match P.iter().reduce(min_distance_to_alpha) {
        Some(x) => {x}
        None => {panic!()} //Idk if this should ever occur
    };

    //corner_2 should be on the other side of alpha, therefor search the other pointset
    let P_for_c2: &Vec<Point>;
    if P1.contains(corner_1) {
        P_for_c2 = P2;
    }
    else if P2.contains(corner_1){
        P_for_c2 = P1;
    }
    else {
        panic!()
    }

    //let min_distance_to_p1 = |p1: &Point, p2: &Point| ( if p1.distance_no_sqrt(p_1) < p2.distance_no_sqrt(p_1) {p1} else {p2} );

    let corner_2 = match P_for_c2.iter().reduce(
        |p1: &Point, p2: &Point| ( if p1.distance_no_sqrt(corner_1) < p2.distance_no_sqrt(corner_1) {p1} else {p2} )
    ) {
        Some(x) => {x}
        None => {panic!()}
    };

    let f = Edge { start: *corner_1, end: *corner_2 };

    //Finds p3 and returns the Simplex
    let t = make_simplex(&f, P);


    let t_inv  = Triangle::new(t.p1, t.p3, t.p2);

    return t_inv
}  

//MakeSimplex selects the point p which minimizes the delauny distance dd (p.335)
//TODO: return option instead if no points in the Halfspace
fn make_simplex(f: &Edge, P: &Vec<Point>) -> Triangle {
  //let min_dd= |p1: &Point, p2: &Point| ( if dd(f, p1) < dd(f, p2) {p1} else {p2} );

  //Only points which are on the outside the already triangulated area
  let mut filtered_P: Vec<Point> = P.iter()
  .filter(|p1| !(*p1).left_side_of_edge(f))
  .cloned()
  .collect(); //TODO: Check !

  //Sort by the distance
  filtered_P.sort_by_key(
    |p1: &Point| ( dd(f, p1) /*.abs()*/ as i32 )
 );

  'outer: for p in &filtered_P {
    let t = Triangle::new(f.start, *p, f.end); 

    for other in &filtered_P {
        if p != other {
            if t.point_inside(*other){
                continue 'outer;
            }
        }
    }

    return t;

  }

  panic!();
}

fn dd(f: &Edge, p: &Point) -> f32 {
  let (center, radius) = p.circumcircle(&f.start, &f.end);
  if in_halfspace(f, p, &center) {
    return radius
  }
  else {
    return -radius
  }
}

fn in_halfspace(f: &Edge, p1: &Point, p2: &Point) -> bool{
    return p1.left_side_of_edge(f) == p2.left_side_of_edge(f)
}

fn pointset_partition(P: &Vec<Point>, alpha: &Edge, P1: &mut Vec<Point>, P2: &mut Vec<Point>) {
    for p in P {
        if p.distance_to_edge(alpha) > 0.0 {
            P1.push(*p);
        }
        else {
            P2.push(*p);
        }
    }
}

fn test_fns(){
    //Test circumcircle
    let p1 = Point::new (0, 2);
    let p2 = Point::new (5, 2);
    let p3 = Point::new (2, 3);
    let p4 = Point::new (9, 2);
    let p5 = Point::new (1, 1);
    let p6 = Point::new (2, 6);
    let p7 = Point::new (2, 0);

    let edge = Edge::new(p1, p2);
    let edge2 = Edge::new(p2, p1);

    println!("Length: {}", p1.length(&p2));
    println!("distance_no_sqrt: {}", p1.distance_no_sqrt(&p2));

    println!("distance_to_edge: {}", p3.distance_to_edge(&edge));
    println!("distance_to_edge: {}", p4.distance_to_edge(&edge));
    println!("distance_to_edge: {}", p5.distance_to_edge(&edge));
    println!("distance_to_edge: {}", p6.distance_to_edge(&edge));

    println!("left: {}", p6.left_side_of_edge(&edge));
    println!("left: {}", p7.left_side_of_edge(&edge));

    let (c, r) = p1.circumcircle(&p2, &p3);
    println!("{} {} {}", c.x, c.y, r);

    println!("{} ", edge == edge2);

    let ddis = dd(&edge, &p3);

    println!("{} ", ddis);
}

fn plot_dd() {
    let imgx = 800;
    let imgy = 600;

    let mut imgbuf = image::RgbImage::new(imgx, imgy);

    let mut points_list: Vec<Point> = Vec::new();
    let edge = Edge::new(
        Point::new(400, 250),
        Point::new(300, 350)
    );


    for x in 0..imgx {
        for y in 0..imgy{
            points_list.push(Point::new(x, y));
        }
    }

    //Only points which are on the outside the already triangulated area
    let mut filtered_P: Vec<Point> = points_list.iter()
        .filter(|p1| !(*p1).left_side_of_edge(&edge))
        .cloned()
        .collect(); //TODO: Check !

    //Sort by the distance
    //filtered_P.sort_by_key(
    //  |p1: &Point| ( dd(edge, p1).abs() as i32 )
    //);




    //for p in points_list{
    for p in filtered_P{
        let dist = dd(&edge, &p);

        let pixel = imgbuf.get_pixel_mut(p.x, p.y);
        if dist > 255.0 {
        *pixel = image::Rgb([255,0,0]);
        }
        else if dist < -255.0 {
            *pixel = image::Rgb([0,0,255]);
        }
        else if dist < 0.0 {
            *pixel = image::Rgb([0,0,(-dist as u8)]);
        }
        else { 
            *pixel = image::Rgb([(dist as u8) ,0,0]);
        }
    }
    imgbuf.save("dd.png").unwrap();

}

fn main() {
    //plot_dd();

    //return;

    let imgx = 800;
    let imgy = 600;

    let mut imgbuf = image::RgbImage::new(imgx, imgy);

    for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
        //let r = (255 * x/imgx) as u8;
        //let g = (255 * y/imgy) as u8;
        
        *pixel = image::Rgb([255,255, 255]);
    }


    let mut points_list: Vec<Point> = Vec::new();
    
    //Add random points

    let amount = 200;
    let mut rng = rand::thread_rng();

    for _i in 0..amount {
        //TODO: Exclude the border region
        let x = rng.gen_range(0..imgx);
        let y = rng.gen_range(0..imgy);

        points_list.push(Point::new(x, y));
    }

    //Plot points
    for p in &points_list {
        let pixel = imgbuf.get_pixel_mut(p.x, p.y);
        *pixel = image::Rgb([0,0,0]);
    }


    //Corner Points
    
    //Sort
    points_list.sort();
    points_list.dedup();

    //Delauny
    let mut afl: Vec<Edge> = Vec::new();
    let edge = Edge { start: Point::new(0, imgy/2), end: Point::new(imgx, imgy/2)};
    let sigma = de_wall(
        points_list,
        Edge::new(Point::new(0, 0), Point::new(0, imgy)),
        Edge::new(Point::new(0, 0), Point::new(imgx, 0)),
        edge,
        &mut afl
    );

    //colourize
    let steps = 500;
        for n in 0..steps {
            let x = edge.start.x as i32 + ((n as f32 / steps as f32) * (edge.end.x as f32 - edge.start.x as f32))as i32;
            let y = edge.start.y as i32 + ((n as f32 / steps as f32) * (edge.end.y as f32 - edge.start.y as f32))as i32;

            let pixel = imgbuf.get_pixel_mut(x as u32, y as u32);
            *pixel = image::Rgb([255,255,0]);
        }

    for triangle in &sigma{

        //Edges in Green
        let steps = 100;
        for n in 0..steps {
            let x = triangle.p1.x as i32 + ((n as f32 / steps as f32) * (triangle.p2.x as f32 - triangle.p1.x as f32))as i32;
            let y = triangle.p1.y as i32 + ((n as f32 / steps as f32) * (triangle.p2.y as f32 - triangle.p1.y as f32))as i32;

            let pixel = imgbuf.get_pixel_mut(x as u32, y as u32);
            *pixel = image::Rgb([0,(255 * n / steps) as u8,0]);
        }
        for n in 0..steps {
            let x = triangle.p2.x as i32 + ((n as f32 / steps as f32) * (triangle.p3.x as f32 - triangle.p2.x as f32))as i32;
            let y = triangle.p2.y as i32 + ((n as f32 / steps as f32) * (triangle.p3.y as f32 - triangle.p2.y as f32))as i32;

            let pixel = imgbuf.get_pixel_mut(x as u32, y as u32);
            *pixel = image::Rgb([0,(255 * n / steps) as u8,0]);
        }
        for n in 0..steps {
            let x = triangle.p3.x as i32 + ((n as f32 / steps as f32) * (triangle.p1.x as f32 - triangle.p3.x as f32))as i32;
            let y = triangle.p3.y as i32 + ((n as f32 / steps as f32) * (triangle.p1.y as f32 - triangle.p3.y as f32))as i32;

            let pixel = imgbuf.get_pixel_mut(x as u32, y as u32);
            *pixel = image::Rgb([0,(255 * n / steps) as u8,0]);
        }


        //Corners in Red
        for point in [triangle.p1, triangle.p2, triangle.p3] {
            let pixel = imgbuf.get_pixel_mut(point.x, point.y);
            *pixel = image::Rgb([255,0,0]);
        }
    }
/*
    if let Some(triangle) = sigma.last(){

        //Edges in Green
        let steps = 100;
        for n in 0..steps {
            let x = triangle.p1.x as i32 + ((n as f32 / steps as f32) * (triangle.p2.x as f32 - triangle.p1.x as f32))as i32;
            let y = triangle.p1.y as i32 + ((n as f32 / steps as f32) * (triangle.p2.y as f32 - triangle.p1.y as f32))as i32;

            let pixel = imgbuf.get_pixel_mut(x as u32, y as u32);
            *pixel = image::Rgb([0,0,255]);
        }
        for n in 0..steps {
            let x = triangle.p2.x as i32 + ((n as f32 / steps as f32) * (triangle.p3.x as f32 - triangle.p2.x as f32))as i32;
            let y = triangle.p2.y as i32 + ((n as f32 / steps as f32) * (triangle.p3.y as f32 - triangle.p2.y as f32))as i32;

            let pixel = imgbuf.get_pixel_mut(x as u32, y as u32);
            *pixel = image::Rgb([0,0,255]);
        }
        for n in 0..steps {
            let x = triangle.p3.x as i32 + ((n as f32 / steps as f32) * (triangle.p1.x as f32 - triangle.p3.x as f32))as i32;
            let y = triangle.p3.y as i32 + ((n as f32 / steps as f32) * (triangle.p1.y as f32 - triangle.p3.y as f32))as i32;

            let pixel = imgbuf.get_pixel_mut(x as u32, y as u32);
            *pixel = image::Rgb([0,0,255]);
        }
    }
    */



    //export
    imgbuf.save("out.png").unwrap();

}

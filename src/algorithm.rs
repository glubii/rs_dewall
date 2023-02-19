// mod geometry;
use crate::geometry::Edge;
use crate::geometry::Point;
use crate::geometry::SpacePartition;
use crate::geometry::Triangle;

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
pub fn de_wall(P: Vec<Point>, sp: &SpacePartition, afl: &mut Vec<Edge>) -> Vec<Triangle> {
    // (d-1) faces intersected by plane alpha
    let mut afl_sigma: Vec<Edge> = Vec::new();

    // (d-1) faces with all of the vertices in P1
    let mut afl_1: Vec<Edge> = Vec::new();

    // (d-1) faces with all of the vertices in P2
    let mut afl_2: Vec<Edge> = Vec::new();

    let mut P1: Vec<Point> = Vec::new();
    let mut P2: Vec<Point> = Vec::new();

    let mut sigma: Vec<Triangle> = Vec::new();

    let mut step_counter: i32 = 0;

    //Divide Points into lower (negative) and upper (positive) groups
    //TODO: Use references instead?
    pointset_partition(&P, &sp.alpha, &mut P1, &mut P2);

    println!("len of P {} {} {}", P.len(), P1.len(), P2.len());

    if P.is_empty() {
        return sigma;
    }

    //Simplex Wall Construction
    if afl.is_empty() {
        let t = make_first_simplex(&P, &sp.alpha, &P1, &P2);

        afl.push(t.e1);
        afl.push(t.e2);
        afl.push(t.e3);
        sigma.push(t);

        println!("First simplex done:");
        println!("1 {} {}", t.p1.x, t.p1.y);
        println!("2 {} {}", t.p2.x, t.p2.y);
        println!("3 {} {}", t.p3.x, t.p3.y);
    }

    for f in afl {
        if P1.contains(&f.start) && P1.contains(&f.end) {
            afl_1.push(*f);
        } else if P2.contains(&f.start) && P2.contains(&f.end) {
            afl_2.push(*f);
        } else {
            //splitting wall intersects triangle
            afl_sigma.push(*f);
        }
    }

    // plot_step(
    //     step_counter,
    //     &sp.alpha,
    //     &P1,
    //     &P2,
    //     &afl_1,
    //     &afl_2,
    //     &afl_sigma,
    //     &sigma,
    // );
    // step_counter += 1;

    //for _face in 0..6 {
    while !afl_sigma.is_empty() {
        // f = Extract(afl_sigma)
        println!("Loop start, {} faces remaining", afl_sigma.len());
        let face_opt = afl_sigma.pop();
        let f = match face_opt {
            Some(x) => x,
            None => {
                panic!()
            } //Idk if this can ever occur
        };

        println!(
            "Looking at face ({} {}), ({} {})",
            f.start.x, f.start.y, f.end.x, f.end.y
        );

        let t = match make_simplex(step_counter, &f, &P) {
            Some(x) => x,
            None => continue,
        };

        // step_counter += 1;

        println!("1 {} {}", t.p1.x, t.p1.y);
        println!("2 {} {}", t.p2.x, t.p2.y);
        println!("3 {} {}", t.p3.x, t.p3.y);

        if true {
            //t != null { TODO: Add option or something to make_simplex
            sigma.push(t);

            for f_apostroph in [t.e1, t.e2, t.e3] {
                if f_apostroph != f {
                    if P1.contains(&f_apostroph.start) && P1.contains(&f_apostroph.end) {
                        //println!("AFL1");
                        //println!("{}", afl_1.len());
                        update(&f_apostroph, &mut afl_1);
                        //println!("{}", afl_1.len());
                    } else if P2.contains(&f_apostroph.start) && P2.contains(&f_apostroph.end) {
                        //println!("AFL2");
                        //println!("{}", afl_2.len());
                        update(&f_apostroph, &mut afl_2);
                        //println!("{}", afl_2.len());
                    } else {
                        //splitting wall intersects triangle
                        //println!("Sigma");
                        //println!("{}", afl_sigma.len());
                        update(&f_apostroph, &mut afl_sigma);
                        //println!("{}", afl_sigma.len());
                    }
                } else {
                    println!("Eq")
                }
                // plot_step(
                //     step_counter,
                //     &sp.alpha,
                //     &P1,
                //     &P2,
                //     &afl_1,
                //     &afl_2,
                //     &afl_sigma,
                //     &sigma,
                // );
                // step_counter += 1;
            }
        }
    }

    //TODO: Why are sp1 and sp2 swapped?
    let (sp1, sp2) = sp.partition();
    if !afl_1.is_empty() {
        sigma.append(&mut de_wall(P1, &sp2, &mut afl_1));
    }
    if !afl_2.is_empty() {
        sigma.append(&mut de_wall(P2, &sp1, &mut afl_2));
    }

    return sigma;
}

//p.336
fn update(f: &Edge, L: &mut Vec<Edge>) {
    if L.contains(f) {
        L.retain(|&x| x != *f);
    } else {
        L.push(*f);
    }
}

//p.335
fn make_first_simplex(P: &Vec<Point>, alpha: &Edge, P1: &Vec<Point>, P2: &Vec<Point>) -> Triangle {
    //Find the point closest to alpha
    //let min_distance_to_alpha = |p1: &Point, p2: &Point| ( if p1.distance_to_edge(alpha) < p2.distance_to_edge(alpha) {p1} else {p2} );

    let corner_1 = match P.iter().reduce(|p1: &Point, p2: &Point| {
        (if p1.distance_to_edge(alpha).abs() < p2.distance_to_edge(alpha).abs() {
            p1
        } else {
            p2
        })
    }) {
        //let p_1= match P.iter().reduce(min_distance_to_alpha) {
        Some(x) => x,
        None => {
            panic!()
        } //Idk if this should ever occur
    };

    //corner_2 should be on the other side of alpha, therefor search the other pointset
    let P_for_c2: &Vec<Point>;
    if P1.contains(corner_1) {
        P_for_c2 = P2;
    } else if P2.contains(corner_1) {
        P_for_c2 = P1;
    } else {
        panic!()
    }

    //let min_distance_to_p1 = |p1: &Point, p2: &Point| ( if p1.distance_no_sqrt(p_1) < p2.distance_no_sqrt(p_1) {p1} else {p2} );

    let corner_2 = match P_for_c2.iter().reduce(|p1: &Point, p2: &Point| {
        (if p1.distance_no_sqrt(corner_1) < p2.distance_no_sqrt(corner_1) {
            p1
        } else {
            p2
        })
    }) {
        Some(x) => x,
        None => {
            panic!()
        }
    };

    let f = Edge {
        start: *corner_1,
        end: *corner_2,
    };

    //Finds p3 and returns the Simplex
    let t = match make_simplex(9999, &f, P) {
        Some(x) => x,
        None => {
            panic!()
        }
    };

    return t;

    // let t_inv  = Triangle::new(t.p1, t.p3, t.p2);

    // return t_inv
}

//MakeSimplex selects the point p which minimizes the delauny distance dd (p.335)
//TODO: return option instead if no points in the Halfspace
fn make_simplex(step_counter: i32, f: &Edge, P: &Vec<Point>) -> Option<Triangle> {
    //let min_dd= |p1: &Point, p2: &Point| ( if dd(f, p1) < dd(f, p2) {p1} else {p2} );

    //Only points which are on the outside the already triangulated area
    let mut filtered_P: Vec<Point> = P
        .iter()
        .filter(|p1| !(*p1).left_side_of_edge(f))
        .cloned()
        .collect(); //TODO: Check !

    //   plot_dd(step_counter, &f, &filtered_P);

    //Sort by the distance
    filtered_P.sort_by_key(|p1: &Point| (dd(f, p1) /*.abs()*/ as i32));

    'outer: for p in &filtered_P {
        if *p == f.start || *p == f.end {
            continue;
        }
        let t = Triangle::new(f.start, *p, f.end);

        println!(
            "Checking point {} {} with delauney distance {}",
            p.x,
            p.y,
            dd(f, p)
        );
        for other in &filtered_P {
            if p != other {
                if t.point_inside(*other) {
                    println!("Point inside best fit, finding next best option");
                    continue 'outer;
                }
            }
        }

        return Some(t);
    }

    return None;
    //   panic!("No solution found");
}

pub fn dd(f: &Edge, p: &Point) -> f32 {
    let (center, radius) = p.circumcircle(&f.start, &f.end);
    if in_halfspace(f, p, &center) {
        return radius;
    } else {
        // return 9999.9
        return -radius;
    }
}

fn in_halfspace(f: &Edge, p1: &Point, p2: &Point) -> bool {
    return p1.left_side_of_edge(f) == p2.left_side_of_edge(f);
}

fn pointset_partition(P: &Vec<Point>, alpha: &Edge, P1: &mut Vec<Point>, P2: &mut Vec<Point>) {
    for p in P {
        if p.left_side_of_edge(alpha) {
            P1.push(*p);
        } else {
            P2.push(*p);
        }
    }
}

fn plot_step(
    step_counter: i32,
    alpha: &Edge,
    P1: &Vec<Point>,
    P2: &Vec<Point>,
    afl_1: &Vec<Edge>,
    afl_2: &Vec<Edge>,
    afl_sigma: &Vec<Edge>,
    sigma: &Vec<Triangle>,
) {
    let imgx = 800;
    let imgy = 600;

    let mut imgbuf = image::RgbImage::new(imgx, imgy);

    //Draw alpha line dividing the sets
    let steps = 500;
    for n in 0..steps {
        let x = alpha.start.x as i32
            + ((n as f32 / steps as f32) * (alpha.end.x as f32 - alpha.start.x as f32)) as i32;
        let y = alpha.start.y as i32
            + ((n as f32 / steps as f32) * (alpha.end.y as f32 - alpha.start.y as f32)) as i32;

        let pixel = imgbuf.get_pixel_mut(x as u32, y as u32);
        *pixel = image::Rgb([255, 255, 0]);
    }

    // Draw points from both sets
    for point in P1 {
        let pixel = imgbuf.get_pixel_mut(point.x, point.y);
        *pixel = image::Rgb([255, 0, 0]);
    }
    for point in P2 {
        let pixel = imgbuf.get_pixel_mut(point.x, point.y);
        *pixel = image::Rgb([255, 255, 0]);
    }

    for triangle in sigma {
        //Edges in Green
        let steps = 100;
        for n in 0..steps {
            let x = triangle.p1.x as i32
                + ((n as f32 / steps as f32) * (triangle.p2.x as f32 - triangle.p1.x as f32))
                    as i32;
            let y = triangle.p1.y as i32
                + ((n as f32 / steps as f32) * (triangle.p2.y as f32 - triangle.p1.y as f32))
                    as i32;

            let pixel = imgbuf.get_pixel_mut(x as u32, y as u32);
            *pixel = image::Rgb([0, (255 * n / steps) as u8, 0]);
        }
        for n in 0..steps {
            let x = triangle.p2.x as i32
                + ((n as f32 / steps as f32) * (triangle.p3.x as f32 - triangle.p2.x as f32))
                    as i32;
            let y = triangle.p2.y as i32
                + ((n as f32 / steps as f32) * (triangle.p3.y as f32 - triangle.p2.y as f32))
                    as i32;

            let pixel = imgbuf.get_pixel_mut(x as u32, y as u32);
            *pixel = image::Rgb([0, (255 * n / steps) as u8, 0]);
        }
        for n in 0..steps {
            let x = triangle.p3.x as i32
                + ((n as f32 / steps as f32) * (triangle.p1.x as f32 - triangle.p3.x as f32))
                    as i32;
            let y = triangle.p3.y as i32
                + ((n as f32 / steps as f32) * (triangle.p1.y as f32 - triangle.p3.y as f32))
                    as i32;

            let pixel = imgbuf.get_pixel_mut(x as u32, y as u32);
            *pixel = image::Rgb([0, (255 * n / steps) as u8, 0]);
        }
    }

    // Active Face lists
    let steps = 50;
    for edge in afl_1 {
        for n in 0..steps {
            let x = edge.start.x as i32
                + ((n as f32 / steps as f32) * (edge.end.x as f32 - edge.start.x as f32)) as i32;
            let y = edge.start.y as i32
                + ((n as f32 / steps as f32) * (edge.end.y as f32 - edge.start.y as f32)) as i32;

            let pixel = imgbuf.get_pixel_mut(x as u32, y as u32);
            *pixel = image::Rgb([127, 0, 127]);
        }
    }
    for edge in afl_2 {
        for n in 0..steps {
            let x = edge.start.x as i32
                + ((n as f32 / steps as f32) * (edge.end.x as f32 - edge.start.x as f32)) as i32;
            let y = edge.start.y as i32
                + ((n as f32 / steps as f32) * (edge.end.y as f32 - edge.start.y as f32)) as i32;

            let pixel = imgbuf.get_pixel_mut(x as u32, y as u32);
            *pixel = image::Rgb([127, 127, 0]);
        }
    }
    for edge in afl_sigma {
        for n in 0..steps {
            let x = edge.start.x as i32
                + ((n as f32 / steps as f32) * (edge.end.x as f32 - edge.start.x as f32)) as i32;
            let y = edge.start.y as i32
                + ((n as f32 / steps as f32) * (edge.end.y as f32 - edge.start.y as f32)) as i32;

            let pixel = imgbuf.get_pixel_mut(x as u32, y as u32);
            *pixel = image::Rgb([0, 0, 127]);
        }
    }

    imgbuf.save(format!("{step_counter}.png")).unwrap();
}

fn plot_dd(step_counter: i32, edge: &Edge, points_list: &Vec<Point>) {
    let imgx = 800;
    let imgy = 600;

    let mut imgbuf = image::RgbImage::new(imgx, imgy);

    for p in points_list {
        let dist = dd(&edge, &p);

        let pixel = imgbuf.get_pixel_mut(p.x, p.y);
        if dist > 255.0 {
            *pixel = image::Rgb([255, 0, 0]);
        } else if dist < -255.0 {
            *pixel = image::Rgb([0, 0, 255]);
        } else if dist < 0.0 {
            *pixel = image::Rgb([0, 255, (-dist as u8)]);
        } else {
            *pixel = image::Rgb([(dist as u8), 255, 0]);
        }
    }
    imgbuf.save(format!("{step_counter}.png")).unwrap();
}

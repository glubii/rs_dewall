use rand::Rng;

mod geometry;
use geometry::Edge;
use geometry::Point;
use geometry::SpacePartition;

mod algorithm;
use algorithm::de_wall;

fn main() {
    let imgx = 1920;
    let imgy = 1080;

    let mut imgbuf = image::RgbImage::new(imgx, imgy);

    for (_x, _y, pixel) in imgbuf.enumerate_pixels_mut() {
        *pixel = image::Rgb([255, 255, 255]);
    }

    let mut points_list: Vec<Point> = Vec::new();

    // Add random points
    let amount = 500;
    let mut rng = rand::thread_rng();

    for _i in 0..amount {
        //TODO: Exclude the border region
        let x = rng.gen_range(0..imgx);
        let y = rng.gen_range(0..imgy);

        let p = Point::new(x, y);
        if !points_list.contains(&p) {
            points_list.push(p);
        }
    }

    //Plot points
    for p in &points_list {
        let pixel = imgbuf.get_pixel_mut(p.x, p.y);
        *pixel = image::Rgb([0, 0, 0]);
    }

    //Corner Points

    //Sort
    points_list.sort();
    points_list.dedup();

    // Apply the Delaunay algorithm
    let mut afl: Vec<Edge> = Vec::new();

    let sp = SpacePartition::new(
        Edge::new(Point::new(0, 0), Point::new(imgx, 0)),
        Edge::new(Point::new(0, 0), Point::new(0, imgy)),
        Edge::new(Point::new(0, imgy / 2), Point::new(imgx, imgy / 2)),
    );

    let sigma = de_wall(points_list, &sp, &mut afl);

    //colourize
    // let steps = 500;
    // for n in 0..steps {
    //     let x = sp.alpha.start.x as i32
    //         + ((n as f32 / steps as f32) * (sp.alpha.end.x as f32 - sp.alpha.start.x as f32))
    //             as i32;
    //     let y = sp.alpha.start.y as i32
    //         + ((n as f32 / steps as f32) * (sp.alpha.end.y as f32 - sp.alpha.start.y as f32))
    //             as i32;

    //     let pixel = imgbuf.get_pixel_mut(x as u32, y as u32);
    //     *pixel = image::Rgb([255, 255, 0]);
    // }

    for triangle in &sigma {
        //Edges in Green
        let steps = 500;
        for n in 0..steps {
            let x = triangle.p1.x as i32
                + ((n as f32 / steps as f32) * (triangle.p2.x as f32 - triangle.p1.x as f32))
                    as i32;
            let y = triangle.p1.y as i32
                + ((n as f32 / steps as f32) * (triangle.p2.y as f32 - triangle.p1.y as f32))
                    as i32;

            let pixel = imgbuf.get_pixel_mut(x as u32, y as u32);
            *pixel = image::Rgb([
                0,
                (255 * x / imgx as i32) as u8,
                (255 * y / imgy as i32) as u8,
            ]);
        }
        for n in 0..steps {
            let x = triangle.p2.x as i32
                + ((n as f32 / steps as f32) * (triangle.p3.x as f32 - triangle.p2.x as f32))
                    as i32;
            let y = triangle.p2.y as i32
                + ((n as f32 / steps as f32) * (triangle.p3.y as f32 - triangle.p2.y as f32))
                    as i32;

            let pixel = imgbuf.get_pixel_mut(x as u32, y as u32);
            *pixel = image::Rgb([
                0,
                (255 * x / imgx as i32) as u8,
                (255 * y / imgy as i32) as u8,
            ]);
        }
        for n in 0..steps {
            let x = triangle.p3.x as i32
                + ((n as f32 / steps as f32) * (triangle.p1.x as f32 - triangle.p3.x as f32))
                    as i32;
            let y = triangle.p3.y as i32
                + ((n as f32 / steps as f32) * (triangle.p1.y as f32 - triangle.p3.y as f32))
                    as i32;

            let pixel = imgbuf.get_pixel_mut(x as u32, y as u32);
            *pixel = image::Rgb([
                0,
                (255 * x / imgx as i32) as u8,
                (255 * y / imgy as i32) as u8,
            ]);
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

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
        // TODO: Exclude the border region
        let x = rng.gen_range(0..imgx);
        let y = rng.gen_range(0..imgy);

        let p = Point::new(x as f32, y as f32);
        if !points_list.contains(&p) {
            points_list.push(p);
        }
    }

    //Plot points
    for p in &points_list {
        let pixel = imgbuf.get_pixel_mut(p.x as u32, p.y as u32);
        *pixel = image::Rgb([0, 0, 0]);
    }

    //Corner Points

    //Sort
    //points_list.sort();
    points_list.dedup();

    // Apply the Delaunay algorithm
    let mut afl: Vec<Edge> = Vec::new();

    let sp = SpacePartition::new(
        Edge::new(Point::new(0.0, 0.0), Point::new(imgx as f32, 0.0)),
        Edge::new(Point::new(0.0, 0.0), Point::new(0.0, imgy as f32)),
        Edge::new(
            Point::new(0.0, imgy as f32 / 2.0),
            Point::new(imgx as f32, imgy as f32 / 2.0),
        ),
    );

    let sigma = de_wall(points_list, &sp, &mut afl);

    for triangle in &sigma {
        // Draw each side of the triangle
        let steps = 500;
        for n in 0..steps {
            let x = triangle.p1.x as i32
                + ((n as f32 / steps as f32) * (triangle.p2.x - triangle.p1.x)) as i32;
            let y = triangle.p1.y as i32
                + ((n as f32 / steps as f32) * (triangle.p2.y - triangle.p1.y)) as i32;

            let pixel = imgbuf.get_pixel_mut(x as u32, y as u32);
            *pixel = image::Rgb([
                0,
                (255 * x / imgx as i32) as u8,
                (255 * y / imgy as i32) as u8,
            ]);
        }
        for n in 0..steps {
            let x = triangle.p2.x as i32
                + ((n as f32 / steps as f32) * (triangle.p3.x - triangle.p2.x)) as i32;
            let y = triangle.p2.y as i32
                + ((n as f32 / steps as f32) * (triangle.p3.y - triangle.p2.y)) as i32;

            let pixel = imgbuf.get_pixel_mut(x as u32, y as u32);
            *pixel = image::Rgb([
                0,
                (255 * x / imgx as i32) as u8,
                (255 * y / imgy as i32) as u8,
            ]);
        }
        for n in 0..steps {
            let x = triangle.p3.x as i32
                + ((n as f32 / steps as f32) * (triangle.p1.x - triangle.p3.x)) as i32;
            let y = triangle.p3.y as i32
                + ((n as f32 / steps as f32) * (triangle.p1.y - triangle.p3.y)) as i32;

            let pixel = imgbuf.get_pixel_mut(x as u32, y as u32);
            *pixel = image::Rgb([
                0,
                (255 * x / imgx as i32) as u8,
                (255 * y / imgy as i32) as u8,
            ]);
        }
    }

    //export
    imgbuf.save("out.png").unwrap();
}

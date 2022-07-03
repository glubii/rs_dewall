use image;
use rand::Rng;

#[derive(PartialEq, Eq, Ord)]
struct Point {
    pub x: u32,
    pub y: u32
}

impl Point {
    fn new (x: u32, y: u32) -> Self {
        Point{x, y}
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

struct Edge;
struct Triangle;

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

    let mut pointsList: Vec<Point> = Vec::new();

    for _i in 0..amount {
        //TODO: Exclude the border region
        let x = rng.gen_range(0..imgx);
        let y = rng.gen_range(0..imgy);

        pointsList.push(Point::new(x, y));

        let pixel = imgbuf.get_pixel_mut(x, y);
        *pixel = image::Rgb([255,255,255]);
    }

    //Corner Points
    
    //Delauny

    //Sort by y value
    pointsList.sort();

    let chunks = pointsList.chunks(3);

    for chunk in chunks {
        let mut it = chunk.iter().peekable();
        while let Some(point) = it.next(){
            if let Some(next_point) = it.peek() {
                println!("{} {}, {} {}", point.x, point.y, next_point.x, next_point.y);

                let res = 100;
                for i in 0..res {
                    let pixel = imgbuf.get_pixel_mut(
                        (point.x as f32 + (i as f32/res as f32) as f32 * (next_point.x as f32 - point.x as f32)) as u32,
                        (point.y as f32 + (i as f32/res as f32) as f32 * (next_point.y as f32 - point.y as f32)) as u32
                    );
                    *pixel = image::Rgb([0,0,0]);

                }
            }
        } 
    }

    //colourize

    //export

    imgbuf.save("out.png").unwrap();
}

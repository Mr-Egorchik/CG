#include <iostream>
#include <png.h>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#define idx(i, j, k, row_len, depth) (i * row_len * depth + j * depth + k)
#define PI 3.14159265

int save_file(std::string filename, int w, int h, int bitdepth, int colortype, unsigned char* picture, int pitch) {
    FILE* fp = fopen(filename.c_str(), "wb");
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    png_infop info_ptr = png_create_info_struct(png_ptr);
    png_set_IHDR(png_ptr,
        info_ptr,
        w,
        h,
        bitdepth,
        colortype,
        PNG_INTERLACE_NONE,
        PNG_COMPRESSION_TYPE_BASE,
        PNG_FILTER_TYPE_BASE);
    png_bytep* row_pointers = (png_bytep*)malloc(sizeof(png_bytep) * h);
    for (int i = 0; i < h; ++i) {
        row_pointers[i] = picture + i * pitch;
    }

    png_init_io(png_ptr, fp);
    png_set_rows(png_ptr, info_ptr, row_pointers);
    png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);

    fclose(fp);
    fp = NULL;
    png_destroy_write_struct(&png_ptr, &info_ptr);
    png_ptr = NULL;
    info_ptr = NULL;
    free(row_pointers);
    row_pointers = NULL;

    return 0;
}

void task_1_1() {
    int w = 256;
    int h = 256;
    unsigned char* img = new unsigned char[w * h];
    for (int i = 0; i < w * h; ++i) {
        img[i] = 0;
    }
    save_file("black.png", w, h, 8, PNG_COLOR_TYPE_GRAY, img, 1 * w);
}

void task_1_2() {
    int w = 256;
    int h = 256;
    unsigned char* img = new unsigned char[w * h];
    for (int i = 0; i < w * h; ++i) {
        img[i] = 255;
    }
    save_file("white.png", w, h, 8, PNG_COLOR_TYPE_GRAY, img, 1 * w);
}

void task_1_3() {
    int w = 256;
    int h = 256;
    unsigned char* img = new unsigned char[w * h * 3];
    for (int i = 0; i < w * h * 3; i += 3) {
        img[i] = 255;
        img[i + 1] = 0;
        img[i + 2] = 0;
    }
    save_file("red.png", w, h, 8, PNG_COLOR_TYPE_RGB, img, 3 * w);
}

void task_1_4() {
    int w = 1024;
    int h = 1024;
    unsigned char* img = new unsigned char[w * h * 3];
    for (int i = 0; i < h; ++i)
        for (int j = 0; j < w; ++j)
            for (int k = 0; k < 3; ++k)
                img[idx(i, j, k, w, 3)] = (i + j + k) % 256;
    save_file("grad.png", w, h, 8, PNG_COLOR_TYPE_RGB, img, 3 * w);
}

struct point {
    double x;
    double y;
    double z;

    point() {
        x = 0;
        y = 0;
        z = 0;
    }

    point(double coords[3]) {
        x = coords[0];
        y = coords[1];
        z = coords[2];
    }
};

struct polygon {
    point x;
    point y;
    point z;

    polygon(point coords[3]) {
        x = coords[0];
        y = coords[1];
        z = coords[2];
    }
};

void read_from_obj(std::string filename, std::vector<point>& res, std::vector<polygon>& polygons) {
    std::ifstream file(filename);
    std::string curr;
    while (getline(file, curr)) {
        if (curr[0] == 'v' && curr[1] == ' ') {
            std::string numb = "";
            double coords[3];
            int pos = 0;
            for (int i = 2; i < curr.size(); ++i) {
                if (curr[i] == ' ') {
                    coords[pos] = std::stod(numb);
                    ++pos;
                    numb = "";
                }
                else {
                    numb.push_back(curr[i]);
                    if (i == curr.size() - 1) {
                        coords[pos] = std::stod(numb);
                    }
                }
            }
            res.push_back(point(coords));
        }
        else if (curr[0] == 'f' && curr[1] == ' ') {
            std::string numb = "";
            point points[3];
            int pos = 0;
            for (int i = 2; i < curr.size(); ++i) {
                if (curr[i] == '/') {
                    points[pos] = res[std::stoi(numb) - 1];
                    ++pos;
                    numb = "";
                    i += numb.length();
                }
                else if (curr[i] == ' ') {
                    numb = "";
                }
                else {
                    numb.push_back(curr[i]);
                }
            }
            polygons.push_back(polygon(points));
        }
    }
}

void line_2_1(int x0, int y0, int x1, int y1, int w, unsigned char*& img) {
    for (float t = 0.0; t < 1.0; t += 0.01) {
        int x = x0 * (1. - t) + x1 * t;
        int y = y0 * (1. - t) + y1 * t;
        for (int k = 0; k < 3; ++k) {
            img[idx(x, y, k, w, 3)] = 255;
        }
    }
}

void task_2_1() {
    int w = 200;
    int h = 200;
    unsigned char* img = new unsigned char[w * h * 3];
    for (int i = 0; i < w * h * 3; ++i) {
        img[i] = 0;
    }
    for (int i = 0; i < 13; ++i) {
        line_2_1(100, 100, 100 + 95 * cos(2 * i * PI / 13), 100 + 95 * sin(2 * i * PI / 13), w, img);
    }

    save_file("task_2_1.png", w, h, 8, PNG_COLOR_TYPE_RGB, img, 3 * w);
}

void line_2_2(int x0, int y0, int x1, int y1, int w, unsigned char*& img) {
    for (int x = x0; x <= x1; x++) {
        float t = (x - x0) / (float)(x1 - x0);
        int y = y0 * (1. - t) + y1 * t;
        for (int k = 0; k < 3; ++k) {
            img[idx(x, y, k, w, 3)] = 255;
        }
    }
}

void task_2_2() {
    int w = 200;
    int h = 200;
    unsigned char* img = new unsigned char[w * h * 3];
    for (int i = 0; i < w * h * 3; ++i) {
        img[i] = 0;
    }
    for (int i = 0; i < 13; ++i) {
        line_2_2(100, 100, 100 + 95 * cos(2 * i * PI / 13), 100 + 95 * sin(2 * i * PI / 13), w, img);
    }

    save_file("task_2_2.png", w, h, 8, PNG_COLOR_TYPE_RGB, img, 3 * w);
}

void line_2_3(int x0, int y0, int x1, int y1, int w, unsigned char*& img) {
    bool steep = false;
    if (std::abs(x0 - x1) < std::abs(y0 - y1)) {
        std::swap(x0, y0);
        std::swap(x1, y1);
        steep = true;
    }
    if (x0 > x1) {
        std::swap(x0, x1);
        std::swap(y0, y1);
    }
    for (int x = x0; x <= x1; x++) {
        float t = (x - x0) / (float)(x1 - x0);
        int y = y0 * (1. - t) + y1 * t;
        if (steep) {
            for (int k = 0; k < 3; ++k) {
                img[idx(y, x, k, w, 3)] = 255;
            }
        }
        else {
            for (int k = 0; k < 3; ++k) {
                img[idx(x, y, k, w, 3)] = 255;
            }
        }
    }
}

void task_2_3() {
    int w = 200;
    int h = 200;
    unsigned char* img = new unsigned char[w * h * 3];
    for (int i = 0; i < w * h * 3; ++i) {
        img[i] = 0;
    }
    for (int i = 0; i < 13; ++i) {
        line_2_3(100, 100, 100 + 95 * cos(2 * i * PI / 13), 100 + 95 * sin(2 * i * PI / 13), w, img);
    }

    save_file("task_2_3.png", w, h, 8, PNG_COLOR_TYPE_RGB, img, 3 * w);
}

void line_2_4(int x0, int y0, int x1, int y1, int w, unsigned char*& img) {
    bool steep = false;
    if (std::abs(x0 - x1) < std::abs(y0 - y1)) {
        std::swap(x0, y0);
        std::swap(x1, y1);
        steep = true;
    }
    if (x0 > x1) { // make it left
        std::swap(x0, x1);
        std::swap(y0, y1);
    }
    int dx = x1 - x0;
    int dy = y1 - y0;
    float derror = std::abs(dy / float(dx));
    float error = 0;
    int y = y0;
    for (int x = x0; x <= x1; x++) {
        if (steep) {
            for (int k = 0; k < 3; ++k) {
                img[idx(y, x, k, w, 3)] = 255;
            }
        }
        else {
            for (int k = 0; k < 3; ++k) {
                img[idx(x, y, k, w, 3)] = 255;
            }
        }
        error += derror;
        if (error > .5) {
            y += (y1 > y0 ? 1 : -1);
            error -= 1;
        }
    }
}

void task_2_4() {
    int w = 200;
    int h = 200;
    unsigned char* img = new unsigned char[w * h * 3];
    for (int i = 0; i < w * h * 3; ++i) {
        img[i] = 0;
    }
    for (int i = 0; i < 13; ++i) {
        line_2_4(100, 100, 100 + 95 * cos(2 * i * PI / 13), 100 + 95 * sin(2 * i * PI / 13), w, img);
    }

    save_file("task_2_4.png", w, h, 8, PNG_COLOR_TYPE_RGB, img, 3 * w);
}

void save_points(unsigned char*& img, int w, int h, std::vector<point>& points, unsigned char point_mask[3]) {
    for (point p : points) {
        for (int k = 0; k < 3; ++k)
            img[idx((int)(-30 * p.y + 500), (int)(30 * p.z + 500), k, w, 3)] = point_mask[k];
    }
    save_file("points.png", w, h, 8, PNG_COLOR_TYPE_RGB, img, 3 * w);
}

void draw_lines(unsigned char*& img, int w, int h, std::vector<polygon>& polygons) {
    for (polygon p : polygons) {
        line_2_4((int)(-30 * p.x.y + 500), (int)(30 * p.x.z + 500), (int)(-30 * p.y.y + 500), (int)(30 * p.y.z + 500), w, img);
        line_2_4((int)(-30 * p.x.y + 500), (int)(30 * p.x.z + 500), (int)(-30 * p.z.y + 500), (int)(30 * p.z.z + 500), w, img);
        line_2_4((int)(-30 * p.z.y + 500), (int)(30 * p.z.z + 500), (int)(-30 * p.y.y + 500), (int)(30 * p.y.z + 500), w, img);
    }
    save_file("lines.png", w, h, 8, PNG_COLOR_TYPE_RGB, img, 3 * w);
}

struct dot {
    double x;
    double y;

    dot() {
        this->x = 0.;
        this->y = 0.;
    }

    dot(double x, double y) {
        this->x = x;
        this->y = y;
    }

    int int_x() {
        return (int)x;
    }

    int int_y() {
        return (int)y;
    }
};

struct triangle {
    dot a;
    dot b;
    dot c;

    triangle(dot a, dot b, dot c) {
        this->a = a;
        this->b = b;
        this->c = c;
    }
};

double* barycentric_coordinates(triangle& tr, dot& d) {
    double* res = new double[3];
    res[0] = ((tr.b.x - tr.c.x) * (d.int_y() - tr.c.y) - (tr.b.y - tr.c.y) * (d.int_x() - tr.c.x)) / ((tr.b.x - tr.c.x) * (tr.a.y - tr.c.y) - (tr.b.y - tr.c.y) * (tr.a.x - tr.c.x));
    res[1] = ((tr.c.x - tr.a.x) * (d.int_y() - tr.a.y) - (tr.c.y - tr.a.y) * (d.int_x() - tr.a.x)) / ((tr.c.x - tr.a.x) * (tr.b.y - tr.a.y) - (tr.c.y - tr.a.y) * (tr.b.x - tr.a.x));
    res[2] = ((tr.a.x - tr.b.x) * (d.int_y() - tr.b.y) - (tr.a.y - tr.b.y) * (d.int_x() - tr.b.x)) / ((tr.a.x - tr.b.x) * (tr.c.y - tr.b.y) - (tr.a.y - tr.b.y) * (tr.c.x - tr.b.x));
    return res;
}

void draw_triangle(triangle& tr, unsigned char*& img, int w, int h, unsigned char point_mask[3]) {
    double xmin = std::min(tr.a.x, std::min(tr.b.x, tr.c.x));
    double xmax = std::max(tr.a.x, std::max(tr.b.x, tr.c.x));
    double ymin = std::min(tr.a.y, std::min(tr.b.y, tr.c.y));
    double ymax = std::max(tr.a.y, std::max(tr.b.y, tr.c.y));
    xmin = xmin < 0 ? 0 : xmin;
    ymin = ymin < 0 ? 0 : ymin;
    xmax = xmax < w ? xmax : w;
    ymax = ymax < h ? ymax : h;
    for (int i = xmin; i < xmax; ++i) {
        for (int j = ymin; j < ymax; ++j) {
            dot d = dot(i, j);
            double* bc = barycentric_coordinates(tr, d);
            if (bc[0] > 0 && bc[1] > 0 && bc[2] > 0) {
                for (int k = 0; k < 3; ++k) {
                    img[idx(i, j, k, w, 3)] = point_mask[k];
                }
            }
            delete(bc);
        }
    }
}

void triangle_task() {
    triangle t = triangle(dot(0., 0.), dot(3., 0.), dot(0., 4.));
    dot d = dot(1, 2);
    int h = 10;
    int w = 10;
    unsigned char* img = new unsigned char[w * h * 3];
    for (int i = 0; i < w * h * 3; ++i) {
        img[i] = 0;
    }
    unsigned char point_mask[3] = { 255, 255, 255 };
    draw_triangle(t, img, w, h, point_mask);
    save_file("triangle.png", w, h, 8, PNG_COLOR_TYPE_RGB, img, 3 * w);
    delete[](img);
    //double* bc = barycentric_coordinates(t, d);
    //std::cout << bc[0] << " + " << bc[1] << " + " << bc[2] << " = " << bc[0] + bc[1] + bc[2] << '\n';
}

void fill_polygons(std::vector<polygon>& polygons, unsigned char*& img, int w, int h) {
    for (polygon p : polygons) {
        unsigned char point_mask[3] = {rand() % 256 , rand() % 256 , rand() % 256 };
        triangle t = triangle(dot(-30 * p.x.y + 500, 30 * p.x.z + 500), dot(-30 * p.y.y + 500, 30 * p.y.z + 500), dot(-30 * p.z.y + 500, 30 * p.z.z + 500));
        draw_triangle(t, img, w, h, point_mask);
    }
    save_file("colored_dog.png", w, h, 8, PNG_COLOR_TYPE_RGB, img, 3 * w);
}

int main()
{
    srand(NULL);
    int h = 1000;
    int w = 1000;
    unsigned char* img = new unsigned char[w * h * 3];
    for (int i = 0; i < w * h * 3; ++i) {
        img[i] = 0;
    }
    std::vector<point> points;
    std::vector<polygon> polygons;
    read_from_obj("dog.obj", points, polygons);
    unsigned char point_mask[3] = { 255, 255, 255 };
    save_points(img, w, h, points, point_mask);
    draw_lines(img, w, h, polygons);
    fill_polygons(polygons, img, w, h);
    delete[](img);
}
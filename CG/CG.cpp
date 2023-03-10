#include <iostream>
#include <png.h>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>
#define idx(i, j, k, row_len, depth) (i * row_len * depth + j * depth + k)
#define PI 3.14159265

using namespace std;

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
    vector<int> polygons = {};
    double norm[3] = {0, 0, 0};
    unsigned char light[3] = { 0, 0, 0 };

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

    point(double x, double y, double z) {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    void add_polygon(int idx) {
        this->polygons.push_back(idx);
    }

    void set_norm(double norm[3]) {
        this->norm[0] = norm[0];
        this->norm[1] = norm[1];
        this->norm[2] = norm[2];
    }

    void set_light() {
        char arg1 = (char)(255 * (norm[2] / sqrt(norm[0] * norm[0] + norm[1] * norm[1] + norm[2] * norm[2]))) % 256;
        light[0] = 0;
        light[1] = 0;
        light[2] = arg1;
    }
};

struct polygon {
    point x;
    point y;
    point z;
    vector<unsigned int> points = {};

    polygon(point coords[3]) {
        x = coords[0];
        y = coords[1];
        z = coords[2];
    }

    polygon(point x, point y, point z) {
        this->x = x;
        this->y = y;
        this->z = z;
    }
};

void read_from_obj(std::string filename, std::vector<point>& res, std::vector<polygon>& polygons) {
    std::ifstream file(filename);
    std::string curr;
    int polygons_counter = 0;
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
            int idx[3];
            int pos = 0;
            for (int i = 2; i < curr.size(); ++i) {
                if (curr[i] == '/') {
                    if (numb != "") {
                        points[pos] = res[std::stoi(numb) - 1];
                        idx[pos] = std::stoi(numb) - 1;
                        ++pos;
                        numb = "";
                        i += numb.length();
                    }
                }
                else if (curr[i] == ' ') {
                    numb = "";
                }
                else {
                    numb.push_back(curr[i]);
                }
            }
            for (int i : idx) {
                res[i].add_polygon(polygons_counter);
            }
            for (auto p : points) {
                p.add_polygon(polygons_counter);
            }
            polygons.push_back(polygon(points));
            for (int i : idx) {
                polygons[polygons.size() - 1].points.push_back(i);
            }
            
            ++polygons_counter;
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
    double norm[3] = {0, 0, 0};

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

    void set_norm(double norm[3]) {
        this->norm[0] = norm[0];
        this->norm[1] = norm[1];
        this->norm[2] = norm[2];
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
        unsigned char point_mask[3] = { rand() % 256 , rand() % 256 , rand() % 256 };
        triangle t = triangle(dot(-30 * p.x.y + 500, 30 * p.x.z + 500), dot(-30 * p.y.y + 500, 30 * p.y.z + 500), dot(-30 * p.z.y + 500, 30 * p.z.z + 500));
        draw_triangle(t, img, w, h, point_mask);
    }
    save_file("colored_dog.png", w, h, 8, PNG_COLOR_TYPE_RGB, img, 3 * w);
}

void find_norm(point& norm, polygon polygon) {
    norm.x = ((polygon.y.y - polygon.x.y) * (polygon.y.z - polygon.z.z) - (polygon.y.z - polygon.x.z) * (polygon.y.y - polygon.z.y));
    norm.y = ((polygon.y.z - polygon.x.z) * (polygon.y.x - polygon.z.x) - (polygon.y.x - polygon.x.x) * (polygon.y.z - polygon.z.z));
    norm.z = ((polygon.y.x - polygon.x.x) * (polygon.y.y - polygon.z.y) - (polygon.y.y - polygon.x.y) * (polygon.y.x - polygon.z.x));
    /*norm.x = ((polygon.y.z - polygon.x.z) * (polygon.y.x - polygon.z.x) - (polygon.y.x - polygon.x.x) * (polygon.y.z - polygon.z.z));
    norm.y = ((polygon.y.x - polygon.x.x) * (polygon.y.y - polygon.z.y) - (polygon.y.y - polygon.x.y) * (polygon.y.x - polygon.z.x));
    norm.z = ((polygon.y.y - polygon.x.y) * (polygon.y.z - polygon.z.z) - (polygon.y.z - polygon.x.z) * (polygon.y.y - polygon.z.y));*/
}

vector<polygon> get_good_polygons(vector<polygon>& polygons) {
    vector<polygon> new_polygons;
    point norm = point();
    for (polygon p : polygons) {
        find_norm(norm, p);
        if (norm.x / sqrt(norm.x * norm.x + norm.y * norm.y + norm.z * norm.z) < 0) {
            new_polygons.push_back(p);
        }
    }
    return new_polygons;
}

void fill_polygons_with_shades(std::vector<polygon>& polygons, unsigned char*& img, int w, int h) {
    point norm = point();
    for (polygon p : polygons) {
        find_norm(norm, p);
        int arg1 = (int)(255 * (norm.x / sqrt(norm.x * norm.x + norm.y * norm.y + norm.z * norm.z))) % 256;
        unsigned char point_mask[3] = { 0, 0, 255 - arg1 };
        triangle t = triangle(dot(-30 * p.x.y + 500, 30 * p.x.z + 500), dot(-30 * p.y.y + 500, 30 * p.y.z + 500), dot(-30 * p.z.y + 500, 30 * p.z.z + 500));
        draw_triangle(t, img, w, h, point_mask);
    }
    save_file("shade_dog.png", w, h, 8, PNG_COLOR_TYPE_RGB, img, 3 * w);
}

void draw_triangle_z_buffer(triangle& tr, polygon& p, double*& z_buffer, unsigned char*& img, int w, int h, unsigned char point_mask[3]) {
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
                double z_0 = bc[0] * p.x.x + bc[1] * p.y.x + bc[2] * p.z.x;
                if (!(z_0 > z_buffer[i * w + j])) {
                    for (int k = 0; k < 3; ++k) {
                        img[idx(i, j, k, w, 3)] = point_mask[k];
                    }
                    z_buffer[i * w + j] = z_0;
                }
            }
            delete[](bc);
        }
    }
}

void fill_polygons_with_z_buffer(std::vector<polygon>& polygons, unsigned char*& img, int w, int h) {
    double* z_buffer = new double[w * h];
    for (int i = 0; i < w * h; ++i) {
        z_buffer[i] = numeric_limits<double>::max();
    }
    point norm = point();
    for (polygon p : polygons) {
        find_norm(norm, p);
        char arg1 = (char)(255 * (norm.x / sqrt(norm.x * norm.x + norm.y * norm.y + norm.z * norm.z))) % 256;
        unsigned char point_mask[3] = { 0, 0, arg1 };
        triangle t = triangle(dot(-30 * p.x.y + 500, 30 * p.x.z + 500), dot(-30 * p.y.y + 500, 30 * p.y.z + 500), dot(-30 * p.z.y + 500, 30 * p.z.z + 500));
        draw_triangle_z_buffer(t, p, z_buffer, img, w, h, point_mask);
    }
    delete[](z_buffer);
    save_file("z_buffered_dog.png", w, h, 8, PNG_COLOR_TYPE_RGB, img, 3 * w);
}

void projective_transform(vector<point>& points, double scale[2], int h, int w, double t[3], vector<point>& new_points) {
    for (auto p : points) {
        point new_point = point(scale[0] * p.x / (p.z + t[2]) + w, scale[1] * p.y / (p.z + t[2]) + h, (p.z + t[2]));
        new_point.polygons = p.polygons;
        new_point.set_norm(p.norm);
        new_point.set_light();
        new_points.push_back(new_point);
    }
}

void upd_save_points(unsigned char*& img, int w, int h, std::vector<point>& points, unsigned char point_mask[3]) {
    for (int i = 0; i < points.size(); ++i) {
        point p = points[i];
        for (int k = 0; k < 3; ++k) {

            img[idx(abs((int)p.x), abs((int)p.y), k, w, 3)] = point_mask[k];
        }
    }
    save_file("upd_points.png", w, h, 8, PNG_COLOR_TYPE_RGB, img, 3 * w);
}

void multiplyMatrices(double mat1[][3], double mat2[][3], double result[][3]) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            result[i][j] = 0;
            for (int k = 0; k < 3; k++) {
                result[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }
}

void rotate(const vector<point>& points, double yaw, double pitch, double roll, vector<point>& output) {
    yaw = yaw * PI / 180.0;
    pitch = pitch * PI / 180.0;
    roll = roll * PI / 180.0;

    // Construct rotation matrices for each axis
    double cosYaw = cos(yaw);
    double sinYaw = sin(yaw);
    double yawMat[3][3] = {
        {cosYaw, 0, sinYaw},
        {0, 1, 0},
        {-sinYaw, 0, cosYaw}
    };

    double cosPitch = cos(pitch);
    double sinPitch = sin(pitch);
    double pitchMat[3][3] = {
        {1, 0, 0},
        {0, cosPitch, sinPitch},
        {0, -sinPitch, cosPitch}
    };

    double cosRoll = cos(roll);
    double sinRoll = sin(roll);
    double rollMat[3][3] = {
        {cosRoll, sinRoll, 0},
        {-sinRoll, cosRoll, 0},
        {0, 0, 1}
    };

    double rotationMatSub[3][3];
    double rotationMat[3][3];
    multiplyMatrices(rollMat, yawMat, rotationMatSub);
    multiplyMatrices(rotationMatSub, pitchMat, rotationMat);

    for (auto p : points) {
        double xRot = rotationMat[0][0] * p.x + rotationMat[0][1] * p.y + rotationMat[0][2] * p.z;
        double yRot = rotationMat[1][0] * p.x + rotationMat[1][1] * p.y + rotationMat[1][2] * p.z;
        double zRot = rotationMat[2][0] * p.x + rotationMat[2][1] * p.y + rotationMat[2][2] * p.z;
        point new_point = point(xRot, yRot, zRot);
        new_point.polygons = p.polygons;
        new_point.set_norm(p.norm);
        new_point.set_light();
        output.push_back(new_point);
    }
}

void set_norm(vector<polygon>& polygons, point& p) {
    double res[3] = { 0, 0, 0 };
    for (int i: p.polygons) {
        polygon polygon = polygons[i];
        res[0] += (((polygon.y.y - polygon.x.y) * (polygon.y.z - polygon.z.z) - (polygon.y.z - polygon.x.z) * (polygon.y.y - polygon.z.y))) / polygons.size();
        res[1] += (((polygon.y.z - polygon.x.z) * (polygon.y.x - polygon.z.x) - (polygon.y.x - polygon.x.x) * (polygon.y.z - polygon.z.z))) / polygons.size();
        res[2] += (((polygon.y.x - polygon.x.x) * (polygon.y.y - polygon.z.y) - (polygon.y.y - polygon.x.y) * (polygon.y.x - polygon.z.x))) / polygons.size();
    }
    p.set_norm(res);
}

void upd_polygons_data(vector<point>& points, vector<polygon>& polygons) {
    for (int i = 0; i < polygons.size(); ++i) {
        polygons[i].x = points[polygons[i].points[0]];
        polygons[i].y = points[polygons[i].points[1]];
        polygons[i].z = points[polygons[i].points[2]];
    }
}

void guro(std::vector<point>& points, vector<polygon>& polygons, unsigned char*& img, int w, int h) {
    double* z_buffer = new double[w * h];
    for (int i = 0; i < w * h; ++i) {
        z_buffer[i] = numeric_limits<double>::max();
    }
    for (polygon p : polygons) {
        point pa = points[p.points[0]];
        point pb = points[p.points[1]];
        point pc = points[p.points[2]];
        dot a = dot(pa.x, pa.y);
        dot b = dot(pb.x, pb.y);
        dot c = dot(pc.x, pc.y);
        triangle tr = triangle(a, b, c);
        a.set_norm(pa.norm);
        b.set_norm(pb.norm);
        c.set_norm(pc.norm);
        double xmin = std::min(tr.a.x, std::min(tr.b.x, tr.c.x));
        double xmax = std::max(tr.a.x, std::max(tr.b.x, tr.c.x));
        double ymin = std::min(tr.a.y, std::min(tr.b.y, tr.c.y));
        double ymax = std::max(tr.a.y, std::max(tr.b.y, tr.c.y));
        xmin = xmin < 0 ? 0 : xmin;
        ymin = ymin < 0 ? 0 : ymin;
        xmax = xmax < w ? xmax : w;
        ymax = ymax < h ? ymax : h;
        double l0 = a.norm[2] / sqrt(a.norm[0] * a.norm[0] + a.norm[1] * a.norm[1] + a.norm[2] * a.norm[2]);
        double l1 = b.norm[2] / sqrt(b.norm[0] * b.norm[0] + b.norm[1] * b.norm[1] + b.norm[2] * b.norm[2]);
        double l2 = c.norm[2] / sqrt(c.norm[0] * c.norm[0] + c.norm[1] * c.norm[1] + c.norm[2] * c.norm[2]);
        for (int i = xmin; i < xmax; ++i) {
            for (int j = ymin; j < ymax; ++j) {
                dot d = dot(i, j);
                double* bc = barycentric_coordinates(tr, d);
                if (bc[0] > 0 && bc[1] > 0 && bc[2] > 0) {
                    double arg1 = bc[0] * l0 + bc[1] * l1 + bc[2] * l2;
                    double z_0 = bc[0] * p.x.z + bc[1] * p.y.z + bc[2] * p.z.z;
                    if(!(z_0 > z_buffer[i * w + j])) {
                        img[idx(i, j, 0, w, 3)] = 0;
                        img[idx(i, j, 1, w, 3)] = 0;
                        img[idx(i, j, 2, w, 3)] = (char)(255 * arg1) % 256;
                        z_buffer[i * w + j] = z_0;
                    }/*
                    printf("arg1 = %f\t color = %d\n", arg1, (char)(255 * arg1) % 255);*/
                }
                delete[](bc);
            }
        }
    }
    save_file("guro_dog.png", w, h, 8, PNG_COLOR_TYPE_RGB, img, 3 * w);
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
    for (int i = 0; i < points.size(); ++i) {
        set_norm(polygons, points[i]);
        points[i].set_light();
    }
    upd_polygons_data(points, polygons);
    unsigned char point_mask[3] = { 255, 255, 255 };
    //save_points(img, w, h, points, point_mask);
    //draw_lines(img, w, h, polygons);
    //fill_polygons(polygons, img, w, h);
    //vector<polygon> new_polys = get_good_polygons(polygons);
    //fill_polygons_with_shades(new_polys, img, w, h);
    //fill_polygons_with_z_buffer(polygons, img, w, h);
    vector<point> new_points;
    double scale[2] = { 2000, 2000 };
    double t[3] = { 0, 0, 50 };
    projective_transform(points, scale, 500, 250, t, new_points);
    vector<point> rotated_new_points;
    rotate(points, 170, -20, -90, rotated_new_points);
    upd_polygons_data(rotated_new_points, polygons);
    for (int i = 0; i < rotated_new_points.size(); ++i) {
        set_norm(polygons, rotated_new_points[i]);
        rotated_new_points[i].set_light();
    }
    upd_polygons_data(rotated_new_points, polygons);
    new_points.clear();
    projective_transform(rotated_new_points, scale, 600, 600, t, new_points);
    //upd_save_points(img, w, h, new_points, point_mask);
    guro(new_points, polygons, img, w, h);
    delete[](img);
}
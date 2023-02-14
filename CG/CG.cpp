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

    point(double coords[3]) {
        x = coords[0];
        y = coords[1];
        z = coords[2];
    }
};

void read_from_obj(std::string filename, std::vector<point>& res) {
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

int main()
{
    task_2_1();
    task_2_2();
    task_2_3();
    task_2_4();
}
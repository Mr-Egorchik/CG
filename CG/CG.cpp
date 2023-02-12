#include <iostream>
#include <png.h>
#define idx(i, j, k, row_len, depth) (i * row_len * depth + j * depth + k)

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

int main()
{
    task_1_1();
    task_1_2();
    task_1_3();
    task_1_4();    
}
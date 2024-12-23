#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>
#include "doslam_class.h"
#include <stdint.h>
#include <sys/time.h>
#include <sys/ioctl.h>

#define REPORT 1

int doslam::matchMask(int x) {

    int result;

    switch (x) {
        case 0:
            result = 0x00e0; // 0000 0000 1110 0000
            break;
        case 7:
            result = 0x0070; // 0000 0000 0111 0000
            break;
        case 6:
            result = 0x0038; // 0000 0000 0011 1000
            break;
        case 5:
            result = 0x001c; // 0000 0000 0001 1100
            break;
        case 4:
            result = 0x000e; // 0000 0000 0000 1110
            break;
        case 3:
            result = 0x0007; // 0000 0000 0000 0111
            break;
        case 2:
            result = 0x8003; // 1000 0000 0000 0011
            break;
        case 1:
            result = 0xc001; // 1100 0000 0000 0001
            break;
        default:
            // Handle default case if needed
            result = 0; // Default value
            break;
    }

    return result;
}

unsigned long doslam::get_time_us() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec * 1000000 + tv.tv_usec;
}

void doslam::readRawImage(const char* filename, unsigned char* buffer, int size) {
    int fd = open(filename, O_RDONLY);
    if (fd < 0) {
        printf("Error opening file\n");
        return;
    }

    if (read(fd, buffer, size) != size) {
        printf("Error reading file\n");
        close(fd);
        return;
    }

    close(fd);
}


void doslam::saveToFile(const char* filename, const void* buf, size_t size) {
    int fd = open(filename, O_RDWR | O_CREAT | O_TRUNC, 0644);  // 以二进制写入模式打开文件
    if (fd == -1) {
        perror("Error opening file");
        return;
    }

    ssize_t bytes_written = write(fd, buf, size);
    if (bytes_written != size) {
        perror("Error writing to file");
        close(fd);
        return;
    }

    close(fd);
    printf("Data saved to file successfully\n");
}

int doslam::init_doslam(){
    // 打开设备文件
    printf("=============1=======================\n");
    printf("fd_ctr: %d\n", fd_ctr);
    printf("fd_data: %d\n", fd_data);
    printf("doslam_image: %p\n", (void*)doslam_image);
    printf("doslam_result: %p\n", (void*)doslam_result);
    fd_ctr = open(DEVICE_CTR_FILE, O_RDWR);
    if (fd_ctr < 0) {
        perror("open");
        return -1;
    }

    // 打开设备文件
    fd_data = open(DEVICE_DATA_FILE, O_RDWR);
    if (fd_data < 0) {
        perror("open");
        return -1;
    }

    doslam_image = (uint8_t *)mmap(NULL, IMAGE_LENGTH, PROT_READ | PROT_WRITE, MAP_SHARED, fd_data, 0);
    if (doslam_image == MAP_FAILED) {
        perror("mmap");
        close(fd_data);
        close(fd_ctr);
        return -1;
    }

    doslam_result = (uint8_t *)mmap(NULL, RSBRIEF_LENGTH, PROT_READ | PROT_WRITE, MAP_SHARED, fd_data, IMAGE_LENGTH);
    if (doslam_result == MAP_FAILED) {
        perror("mmap");
        close(fd_data);
        close(fd_ctr);
        return -1;
    }
    printf("--------------Initial successfully--------------\n");
    printf("=============2=======================\n");
    printf("fd_ctr: %d\n", fd_ctr);
    printf("fd_data: %d\n", fd_data);
    printf("doslam_image: %p\n", (void*)doslam_image);
    printf("doslam_result: %p\n", (void*)doslam_result);
    return 0;
}

void doslam::init_user_data(doslam_ioctl_data *user_data){
    printf("Entering init_user_data\n");
    printf("user_data: %p\n", (void*)user_data);
    printf("siezof(*user_data) = %d \n", sizeof(*user_data));
    if (user_data == nullptr) {
        fprintf(stderr, "Error: user_data is null!\n");
        return; // 或者处理错误
    }
    user_data->sizeInRow = 480;
    user_data->sizeInCol = 640;
    user_data->threshold = 30;
    user_data->topNum = 256;
    user_data->thresholdInit = -1;
    user_data->dmaImageReadAddr = 0;
    user_data->dmaImageWriteAddr = 0x4B400;
    user_data->dmaOrbWriteAddr = IMAGE_LENGTH;
}

void doslam::next_user_data(doslam_ioctl_data *user_data){
    user_data->sizeInRow = user_data->sizeOutRow;
    user_data->sizeInCol = user_data->sizeOutCol;
    user_data->topNum = user_data->topNum * 4 / 5;
    uint32_t temp = user_data->dmaImageReadAddr;
    user_data->dmaImageReadAddr = user_data->dmaImageWriteAddr;
    user_data->dmaImageWriteAddr = temp;

    user_data->dmaOrbWriteAddr = user_data->dmaOrbWriteAddr + user_data->outputLength * 64;
}

void doslam::report(doslam_ioctl_data *user_data){
    int res = ioctl(fd_ctr, 2222, &user_data);

    printf("res %d\n",res);

    printf("dmaImageReadReady %d\n", user_data->dmaImageReadReady);
    printf("dmaImageWriteReady %d\n",user_data->dmaImageWriteReady);
    printf("dmaOrbWriteReady %d\n",  user_data->dmaOrbWriteReady);
    printf("inputLength %d\n",       user_data->inputLength);
    printf("outputLength %d\n",      user_data->outputLength);
    printf("sizeOutRow %d\n",        user_data->sizeOutRow);
    printf("sizeOutCol %d\n",        user_data->sizeOutCol);
}

void doslam::work(doslam_ioctl_data *user_data){
    int res = ioctl(fd_ctr, DOSLAM_IOCTL_CMD_STAR_WORK, user_data);

    int k = 0;
    while(ioctl(fd_ctr, DOSLAM_IOCTL_CMD_READ_STATE, user_data) != 1){
        k++;
        // usleep(1);
        if(k > 1000){
            printf("ka le\n");
            break;
        }
    }
    #ifdef OURPUT_FIRE
        sprintf(str, "/home/ubuntu/doslam/img/1_r%d_c%d.raw", user_data->sizeOutRow, (user_data->sizeOutCol+7)/8*8);
        saveToFile(str, (void *)(doslam_image + user_data->dmaImageWriteAddr), user_data->sizeOutRow * (((user_data->sizeOutCol+7)/8)*8));
        sprintf(str, "/home/ubuntu/doslam/img/1_r%d_c%d.raw.data", user_data->sizeInRow, (user_data->sizeInCol+7)/8*8);
        saveToFile(str, (void *)(doslam_result + user_data->dmaOrbWriteAddr - IMAGE_LENGTH), user_data->outputLength * 64);
    #endif

    #ifdef REPORT
        report(user_data);
    #endif
}


doslam::doslam(/* args */) : fd_ctr(-1), fd_data(-1), doslam_image(nullptr), doslam_result(nullptr)
{
    init_doslam();
    for(int i = 0; i < 4; ++i) {
        orb_addr[i] = nullptr;
        orb_unm[i] = 0;
    }
}

doslam::~doslam()
{
}
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>
#include "doslam.h"
#include <stdint.h>

//================new========================================
#define printfl(x,y) printf("%x,%x\n",x,y);

class doslam
{
private:
    /* data */
public:
    uint8_t *doslam_image;
    uint8_t *doslam_result;

    char str[100];
    uint8_t * orb_addr[4]; // 存储4层orb特征的内存写地址
    uint orb_unm[4]; // 记录每一层输出的特征点个数，每个64B，存储格式为{0，32B->des, 4B->Y, 4B->X, 4B->score}
    int fd_ctr;
    int fd_data;
    int matchMask(int x);
    unsigned long get_time_us();
    void readRawImage(const char* filename, unsigned char* buffer, int size);
    void saveToFile(const char* filename, const void* buf, size_t size);
    int init_doslam();
    void init_user_data(doslam_ioctl_data *user_data);
    void next_user_data(doslam_ioctl_data *user_data);
    void work(doslam_ioctl_data *user_data);
    void report(doslam_ioctl_data *user_data);


    doslam_ioctl_data user_data;
    doslam(/* args */);
    ~doslam();
};


#include <iostream>
#include "doslam_class.h"
#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>
#include "ORBextractor.h"
using namespace std;

#if USE_ORBSLAM2
    doslam kr260;
#endif


int main(){
    cv::Mat im;
    im = cv::imread("/home/ubuntu/TUM_dataset/rgbd_dataset_freiburg1_xyz/rgb/1305031102.175304.png");
    cv::Mat mImGray = im;
    int mbRGB = 1;
   if(mImGray.channels()==3)
    {
        if(mbRGB)
            cv::cvtColor(mImGray,mImGray,CV_RGB2GRAY);
        else
            cv::cvtColor(mImGray,mImGray,CV_BGR2GRAY);
    }
    else if(mImGray.channels()==4)
    {
        if(mbRGB)
            cv::cvtColor(mImGray,mImGray,CV_RGBA2GRAY);
        else
            cv::cvtColor(mImGray,mImGray,CV_BGRA2GRAY);
    }

    std::vector<cv::KeyPoint> mvKeys;
    cv::Mat mDescriptors;
    ORB_SLAM2::ORBextractor* mpORBextractorLeft = new ORB_SLAM2::ORBextractor(1000, 1.25, 4, 20, 7);;
    (*mpORBextractorLeft)(mImGray,cv::Mat(),mvKeys,mDescriptors);
    return 0;
}

#include<algorithm>

#ifndef OPENCV_H
#define OPENCV_H
	#include <opencv2/opencv.hpp>
	using namespace cv;
#endif
	using namespace std;

class Reference {
public:
	int getPerf(int idx);
	int clearPerf(int idx);

	void Resize_Linear(Mat &matSrc, Mat &matDst);
	void Resize_Linear_Test(Mat &matSrc, Mat &matDst);
	void Resize_Linear_Test_1_2(Mat &matSrc, Mat &matDst);
	vector<int> Windows_Test(Mat &matSrc, int size_h, int size_w, int data_num);
	vector<Mat> ComputePyramid(cv::Mat image);

	void RSBrief(const Mat& img, const KeyPoint& kpt, const Point* pattern, uchar* desc);
	void FAST_SCORE(InputArray _img, std::vector<int>& keypoints, int threshold, std::vector<int>& temp);
	void FAST_Test(InputArray _img, std::vector<KeyPoint>& keypoints, int threshold, std::vector<int>& temp);
	void GaussianBlur(Mat &image, Mat &my_res);
//	void ReflectionFillWindow(Mat &_img, std::vector<int64_t>& keypoints);
	void ReflectionFillWindow(Mat &_img, std::vector<int64_t>& keypoints,int h, int w);
	void GaussianBlur(Mat &image, Mat &cv_res, Mat &my_res);
	void FastData(Mat &image, std::vector<int>& keypoints);
	void FAST_Detection(InputArray _img, std::vector<int>& keypoints, int threshold, std::vector<int>& temp);
	void CornerScore(InputArray _img, std::vector<int>& keypoints, int threshold, std::vector<int>& temp);
	void NMS(InputArray _img, std::vector<int>& keypoints, int threshold, std::vector<int>& temp);
	void NMS_cut(InputArray _img, std::vector<int>& keypoints, int threshold, std::vector<int>& temp);
	void noPaddingWindow(Mat &_img, std::vector<int>& keypoints);
	void IC_Angle(Mat &image, Point2f pt, int* mx, int* my);
	void BRIEF(Mat& image, Point2f pt, std::vector<int>& res);
	void tan2(int y, int x, int& res);
	void tan2_q(int y, int x, int& res, int acc);
	void rotate(float angle, int idx);
	void RSBrief(Mat& img, std::vector<KeyPoint>& keypoints, vector<vector<uint32_t>> &descriptors);
};

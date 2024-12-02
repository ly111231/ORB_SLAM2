#ifndef OPENCV_H
#define OPENCV_H
#include <opencv2/opencv.hpp>

using namespace cv;
#endif

#include <fstream>
using namespace std;

const String DIRPATH = "/home/qingyu/orbslam/data/test/";

class Verification {
public:
	void Output_Data_File(Mat &image, String fileName);
	void Output_Data_File(Mat &image, String fileName,int colNum);
	void Output_Data_File(vector<int> data, String fileName, int colNum);
	void Output_Data_File(vector<int> data, String fileName, int colNum, int Size);
	void Output_Data_File8(vector<int64_t> data, String fileName);
	void output_FeaturePoint_File(vector<KeyPoint> &keypoints, String fileName);
	void Resize_Linear_Verification(Mat &image);
	void Windows_Test_Verification(Mat &matSrc, int size_h, int size_w, int data_num);
	void BitwidthConversion_Test_Verification(Mat &matSrc);
	void FAST_Test_Verification(Mat &matSrc);
	void ComputePyramid(Mat &image);

    void Top_orb(Mat& image);

	unordered_map<int, int> subtractArrays(int arr1[], int arr2[], int n);
	void subtractArrays(const vector<int>& arr1, const vector<int>& arr2);

    vector<int> matToVecFill(Mat &mat);
	vector<int> matToVec(Mat &mat);

	void TestReflectionFillWindow(Mat &image,int h,int w);
	void ReflectionGaussianBlur(Mat &image);
	void LoadImages(const string &strFile, vector<string> &vstrImageFilenames, vector<double> &vTimestamps);
	void mul_Image_FAST_Test_Verification();
	void FAST_Detection_Test(InputArray _img, int threshold);
	void CornerScore_Test(InputArray _img, int threshold);
	void NMS_Test(InputArray _img, int threshold);
	void FAST_Test(Mat _img, int threshold);
	void mul_Image_FAST_Test_Time_Verification(int time);
	void DataGenerateIC_Angle(Mat& image, Point2f pt, std::vector<int>& DataInput);
	void IC_Angle_Test(Mat img, int n);
	void BRIEF_Test(Mat& img, int n);
	void tan2_Test(int n);
	void rotate_Test(int n);
	void rotate_Test1(int n);
	void RSBrief_Test(Mat& image, std::vector<KeyPoint>& keypoints, vector<vector<uint32_t>> &descriptors);
};

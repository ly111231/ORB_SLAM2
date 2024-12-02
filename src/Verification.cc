#include "Verification.h"
#include "Reference.h"
#include <chrono>
#include <thread>
#include <iostream>
//这个代码的目标是为了验证，参考代码Resize_Linear和框架代码resize得到的结果能否相同，相差多少。
void Verification::Resize_Linear_Verification(Mat &image) {
	//打印输出的目录
	String ResizePath = "Resize/";
	//定义参考模型的原始数据的输
	int dcol = (image.cols / 5) * 4 + ((image.cols % 5) == 0 ? 0 : (image.cols % 5) - 1);
    int drow = (image.rows / 5) * 4 + ((image.rows % 5) == 0 ? 0 : (image.rows % 5) - 1);
    printf("dest col %d, row %d\n",dcol,drow);
    cv::Size dsize((int)(dcol), (int)(drow));

	Mat matSrcR = image.clone();
	Mat matDstR = Mat::zeros(dsize, image.type());
	Reference ref;
	//使用参考模型得到结果
	ref.Resize_Linear_Test_1_2(matSrcR, matDstR);
	//定义框架方法的原始数据的输入
	Mat matSrcD = image.clone();
	Mat matDstD;
	resize(matSrcD, matDstD, dsize, 0.8, 0.8, INTER_LINEAR);
	//定义resize方法
	//打印结果
	vector<int> vecDstD = matToVec(matDstD);
	vector<int> vecDstR = matToVec(matDstR);
	imshow("参考输出", matDstR);

	subtractArrays(vecDstD, vecDstR);
	Output_Data_File(matToVecFill(image), ResizePath + "SourceDataIn.txt", 8);
	Output_Data_File(matToVecFill(matDstR), ResizePath + "ReferenceDataOut.txt", 8);
	Output_Data_File(matToVecFill(matDstD), ResizePath + "OpencvDataOut.txt", 8);
}


void Verification::Windows_Test_Verification(Mat &matSrc, int size_h, int size_w, int data_num) {
	String ResizePath = "Windows\\";
	// 定义裁剪区域（左上角坐标为 (x,y)，宽度为 w，高度为 h）
	//int x = 100, y = 100, w = 200, h = 200;
	//cv::Rect cropRegion(x, y, w, h);
	//Mat matSrcR = matSrc(cropRegion);
	Mat matSrcR = matSrc.clone();
	Reference ref;
	//使用参考模型得到结果
	vector<int> vecSrcR = matToVec(matSrcR);
	vector<int> vecDstR = ref.Windows_Test(matSrcR, size_h, size_w, data_num);
	Output_Data_File(vecSrcR, ResizePath + "ReferenceDataIn.txt", data_num);
	Output_Data_File(vecDstR, ResizePath + "ReferenceDataOut.txt", data_num);
}

void Verification::BitwidthConversion_Test_Verification(Mat &matSrc) {
	String BitwidthConversionPath = "BitwidthConversion\\";
	vector<int> data = matToVec(matSrc);
	cout << data.size() << endl;
	cout << data.size()%8 << endl;
	Output_Data_File(data, BitwidthConversionPath + "referenceDataIn.txt", 8);
	Output_Data_File(data, BitwidthConversionPath + "referenceDataOut.txt", 10);
}

void Verification::FAST_Test_Verification(Mat &matSrc) {
	 
	String FastPath = "Fast/";
	
	vector<KeyPoint> fast_opencv, fast_reference;
	Reference ref;
	int fastInit = 20;
	FAST(matSrc, fast_opencv, fastInit, true);
	std::vector<int> myVector(9);
	ref.FAST_Test(matSrc, fast_reference, fastInit, myVector);
	cout << "fast:opencv" << endl;
	ofstream ofs1(DIRPATH + FastPath + "OpencvDataOut.txt");
	for (auto kp : fast_opencv) {
		ofs1 << "x:" << kp.pt.x << " y:" << kp.pt.y << " sorce:" << kp.response << endl;
	}
	Mat disimg_opencv;
	drawKeypoints(matSrc, fast_opencv, disimg_opencv);
	imshow("a", disimg_opencv);
	Mat disimg;
	drawKeypoints(matSrc, fast_reference, disimg);
	imshow("b", disimg);
	//waitKey();
	ofstream ofs2(DIRPATH + FastPath + "ReferenceDataOut.txt");
	cout << "fast:reference" << endl;
	for (auto kp : fast_reference) {
		ofs2 << "x:" << kp.pt.x << " y:" << kp.pt.y << " sorce:" << kp.response << endl;
	}
//	std::vector<int64_t> keypoints;
//	ref.ReflectionFillWindow(matSrc, keypoints);
//	Output_Data_File8(keypoints, FastPath + "ReferenceDataIn.txt");
    Output_Data_File(matToVecFill(matSrc), FastPath + "SourceDataIn.txt", 8);
    output_FeaturePoint_File(fast_reference, FastPath + "opencvDataOutFp.coe");
	int errorCount = 0;
}



//将mat类型的图片数据输出到文件中。
void Verification::Output_Data_File(Mat &image, String fileName) {
    Output_Data_File(image, fileName, 8);
}

void Verification::Output_Data_File(Mat &image, String fileName,int colNum) {
    ofstream ofs(DIRPATH + fileName);
    int w = image.cols;
    int h = image.rows;
    int dims = image.channels();
    for (int row = 0; row < h; row++) {
        for (int col = 0; col < w; col++) {
            for (int dim = 0; dim < dims; dim++) {
                ofs << hex << setfill('0') << setw(2) << (int)image.at<Vec3b>(row, col)[dim];
            }
            if (row == h - 2 && col == w - 2)
                continue;
            if(col % colNum == 7){
                ofs << endl;
            }
        }
        int a = colNum - (w % colNum);
        if(w % colNum != 0){
            while(a--){
                ofs << hex << setfill('0') << setw(2) << 0;
            }
            ofs << endl;
        }
    }
}

void Verification::Output_Data_File(vector<int> data, String fileName, int colNum) {
	Output_Data_File(data, fileName, colNum, 2);
}

void Verification::output_FeaturePoint_File(vector<KeyPoint> &keypoints, String fileName){
	ofstream ofs(DIRPATH + fileName);
	for (auto kp : keypoints) {
		ofs << hex << setfill('0') << setw(4) << (int)kp.pt.x;
		ofs << hex << setfill('0') << setw(4) << (int)kp.pt.y;
		ofs << hex << setfill('0') << setw(2) << (int)kp.response;
		ofs << endl;
    }
}

void Verification::Output_Data_File8(vector<int64_t> data, String fileName){
	string file_name = DIRPATH + fileName;
	ofstream ofs(file_name);
	for (int i = 0; i < data.size(); i++) {
		ofs << hex << setfill('0') << setw(16) << data[i] << endl;
	}
}
//将vector类型的数据输出到文件中。
void Verification::Output_Data_File(vector<int> data, String fileName, int colNum, int Size) {
	string file_name = DIRPATH + fileName;

	cout << "1 row data size is:" << colNum << endl;
	ofstream ofs(file_name);
	int temp[20] = { 0 };
//	if (data.size % colNum != 0) {
//		cout << "debug Output_Data_File" << endl;
//	}
	cout << " output data to " << DIRPATH + fileName << endl;
	while (data.size() % colNum != 0) {
		cout << "padding zero" << endl;
		data.push_back(0);
	}
	for (int i = 0; i < data.size(); i += colNum) {
		for (int a = colNum - 1; a > 0; a--) {
			ofs << hex << setfill('0') << setw(Size) << data[i + a];
		}

		if (i != data.size() - 1) {
			ofs << hex << setfill('0') << setw(Size) << data[i] << endl;
		}
		else {
			ofs << hex << setfill('0') << setw(Size) << data[i];
		}
	}
	cout << "data output done rows:" << data.size() / colNum << endl;
}

//对两个int数组输入进行对比，输出两个文件相差的数值，并统计，按统计大小排序
unordered_map<int, int> Verification::subtractArrays(int arr1[], int arr2[], int n) {
	unordered_map<int, int> count;
	cout << "Result: ";
	for (int i = 0; i < n; i++) {
		int result = arr1[i] - arr2[i];
		cout << result << " ";
		count[result]++;
	}
	cout << endl;
	vector<pair<int, int>> counts(count.begin(), count.end());
	sort(counts.begin(), counts.end(), [](pair<int, int>& a, pair<int, int>& b) {
		return a.second > b.second;
	});
	for (auto it : counts) {
		if (it.second > 1) {
			cout << "Result " << it.first << " occurs " << it.second << " times." << endl;
		}
	}
	return count;
}

//对两个int数组输入进行对比，输出两个文件相差的数值，并统计，按统计大小排序
void Verification::subtractArrays(const vector<int>& arr1, const vector<int>& arr2) {
	unordered_map<int, int> count;
	//判断大小
	if (arr1.size() != arr2.size()) {
		cout << "Unequal size arr1 is " << arr1.size() << " arr2 is " << arr2.size() << endl;
	}
	else {
		cout << "Equal size arr1 and arr2 is " << arr1.size() << endl;
	}
	size_t n = min(arr1.size(), arr2.size());
	//计算差值
	for (int i = 0; i < n; i++) {
		int result = arr1[i] - arr2[i];
		count[result]++;
	}
	//按插值总量排序
	vector<pair<int, int>> counts(count.begin(), count.end());
	sort(counts.begin(), counts.end(), [](pair<int, int>& a, pair<int, int>& b) {
		return a.second > b.second;
	});
	//统计正确结果
	cout << "Number of correct results is " << count[0] << endl;
	cout << "The total number is" << n << endl;
	cout << "The accuracy rate is" << (double)count[0] / n << endl;
	//统计错误结果
	for (auto it : counts) {
		if (it.second > 1 && it.first != 0) {
			cout << "Error results " << it.first << " occurs " << it.second << " times.Error ratio is" << (double)it.second / n << endl;
		}
	}
}




//对mat的数值进行类型转换变为int数组,不进行填充
vector<int> Verification::matToVec(Mat &mat) {
	int channels = mat.channels();
	uchar* dataDst = mat.data;
	vector<int> arr1;
	for (int i = 0; i < mat.rows; i++) {
		// 获取第 i 行的指针
		uchar* ptr = mat.ptr<uchar>(i);
		for (int j = 0; j < mat.cols; j++) {
			// 获取第 j 列的像素值
			for (int k = 0; k < mat.channels(); k++) {
				arr1.push_back((int)ptr[j * mat.channels() + k]);
			}
		}
	}
	return arr1;
}

//对mat的数值进行类型转换变为int数组,如果不是8的倍数，那么填充0,返回一个Vector数组，
vector<int> Verification::matToVecFill(Mat &mat) {
    int channels = mat.channels();
    uchar* dataDst = mat.data;
    vector<int> arr1;
    for (int i = 0; i < mat.rows; i++) {
        // 获取第 i 行的指针
        uchar* ptr = mat.ptr<uchar>(i);
        for (int j = 0; j < mat.cols; j++) {
            // 获取第 j 列的像素值
            for (int k = 0; k < mat.channels(); k++) {
                arr1.push_back((int)ptr[j * mat.channels() + k]);
            }
        }
        while(arr1.size() % 8 != 0){
            arr1.push_back(0);
        }
    }
    return arr1;
}


void Verification::ComputePyramid(Mat &image)
{
	int EDGE_THRESHOLD = 18;
	Size sz(cvRound((float)image.cols), cvRound((float)image.rows));
	Size wholeSize(sz.width + EDGE_THRESHOLD * 2, sz.height + EDGE_THRESHOLD * 2);
	Mat temp(wholeSize, image.type()), masktemp;
	Mat mvImagePyramid = temp(Rect(EDGE_THRESHOLD, EDGE_THRESHOLD, sz.width, sz.height));

	copyMakeBorder(image, temp, EDGE_THRESHOLD, EDGE_THRESHOLD, EDGE_THRESHOLD, EDGE_THRESHOLD,
		BORDER_REFLECT_101);

	imshow("原图", image);
	imshow("padding图", temp);
	imshow("特征金子塔图像", mvImagePyramid);
	imshow("fast图像", mvImagePyramid.rowRange(EDGE_THRESHOLD, EDGE_THRESHOLD + 400).colRange(EDGE_THRESHOLD, EDGE_THRESHOLD + 400));
}


void Verification::TestReflectionFillWindow(Mat &image,int h,int w) {
	if (image.cols % 8 != 0) {
		cout << "ReflectionFillWindow error" << endl;
		return;
	}
	Reference ref;
	std::vector<int64_t> vecDstD;
	ref.ReflectionFillWindow(image, vecDstD,h ,w);
	String ReflectionFillPath = "ReflectionFillWindow/";
	vector<int> vecDstR = matToVec(image);
	Output_Data_File(vecDstR, ReflectionFillPath + "ReferenceDataIn.txt", 8);
	Output_Data_File8(vecDstD, ReflectionFillPath + "ReferenceDataOut.txt");
}


void Verification::ReflectionGaussianBlur(Mat &image) {
//	if (image.cols % 8 != 0) {
//		cout << "ReflectionGaussianBlur error" << endl;
//		return;
//	}
	Reference ref;
	//获得输入格式的数据
	std::vector<int> vecDstR;
//	ref.ReflectionFillWindow(image, vecDstR);
	//std::vector<int> vecDstR;
	Mat cvImage, myImage;
	ref.GaussianBlur(image, cvImage, myImage);
	String ReflectionFillPath = "GaussianBlur/";
//	vector<int> vecDstR = matToVec(image);
	std::vector<int> vecDstDcv = matToVecFill(cvImage);
	std::vector<int> vecDstDmy = matToVecFill(myImage);
    Output_Data_File(matToVecFill(image), ReflectionFillPath + "SourceDataIn.txt", 8);
	Output_Data_File(vecDstR, ReflectionFillPath + "ReferenceDataIn.txt", 8);
	Output_Data_File(vecDstDcv, ReflectionFillPath + "ReferenceDataOut.txt", 8);

	subtractArrays(vecDstDcv, vecDstDmy);
}


void Verification::LoadImages(const string &strFile, vector<string> &vstrImageFilenames, vector<double> &vTimestamps)
{
	ifstream f;
	f.open(strFile.c_str());
	cout << strFile << endl;
	// skip first three lines
	string s0;
	getline(f, s0);
	getline(f, s0);
	getline(f, s0);
	cout << "11" << endl;
	cout << s0 << endl;
	while (!f.eof())
	{
		string s;
		getline(f, s);
		if (!s.empty())
		{
			stringstream ss;
			ss << s;
			double t;
			string sRGB;
			ss >> t;
			vTimestamps.push_back(t);
			ss >> sRGB;
			vstrImageFilenames.push_back(sRGB);
		}
	}
}

void Verification::mul_Image_FAST_Test_Verification() {
	vector<string> vstrImageFilenames;
	vector<double> vTimestamps;
	Verification ver;
	Reference ref;
	cout << "11" << endl;
	String FastPath = "F:/Data/orbSlamData/rgbd_dataset_freiburg1_xyz";
	ver.LoadImages(FastPath + "/rgb.txt", vstrImageFilenames, vTimestamps);
	cout << "11" << endl;
	Mat im;
	std::vector<int> myVector(9);
	for (int ni = 0; ni < vstrImageFilenames.size(); ni++) {
		im = cv::imread(FastPath + "/" + vstrImageFilenames[ni], IMREAD_GRAYSCALE);
		int fastInit = 7;
		vector<KeyPoint> fast_opencv, fast_reference;
		//FAST(matSrc, fast_opencv, fastInit, true);
		ref.FAST_Test(im, fast_reference, fastInit, myVector);
		Mat disimg_opencv;
		drawKeypoints(im, fast_reference, disimg_opencv);
		cv::imshow("Image", disimg_opencv);

		// 等待一段时间
		cv::waitKey(1);
	}
	for (int i = 0; i < 9; i++) {
		cout << "8 col sum:" << i << "is count:" << myVector[i] <<" data" << myVector[i] / vstrImageFilenames.size() << endl;
	}
//	cout << "fast:opencv" << endl;
//	ofstream ofs1(DIRPATH + FastPath + "OpencvDataOut.txt");
//	for (auto kp : fast_opencv) {
//		ofs1 << "x:" << kp.pt.x << " y:" << kp.pt.y << " sorce:" << kp.response << endl;
//	}
//	Mat disimg_opencv;
//	drawKeypoints(matSrc, fast_opencv, disimg_opencv);
//	imshow("a", disimg_opencv);
//	Mat disimg;
//	drawKeypoints(matSrc, fast_reference, disimg);
//	imshow("b", disimg);
//	waitKey();
//	ofstream ofs2(DIRPATH + FastPath + "referenceDataOut.txt");
//	cout << "fast:reference" << endl;
//	for (auto kp : fast_reference) {
//		ofs2 << "x:" << kp.pt.x << " y:" << kp.pt.y << " sorce:" << kp.response << endl;
//	}
//	int errorCount = 0;
}

void Verification::mul_Image_FAST_Test_Time_Verification(int time) {
	vector<string> vstrImageFilenames;
	vector<double> vTimestamps;
	int fastInit = 7;
	Verification ver;
	Reference ref;
	String FastPath = "F:/Data/orbSlamData/rgbd_dataset_freiburg1_xyz";
	ver.LoadImages(FastPath + "/rgb.txt", vstrImageFilenames, vTimestamps);
	Mat im;
	std::vector<int> myVector(9);
	int count = 0;
	
	for (int ni = 0; ni < vstrImageFilenames.size(); ni++) {
		im = cv::imread(FastPath + "/" + vstrImageFilenames[ni], IMREAD_GRAYSCALE);
		
		if (ni == 0) {
			cout << im.rows<<" "<<im.cols << endl;
		}
		vector<KeyPoint> fast_opencv, fast_reference;
		//FAST(matSrc, fast_opencv, fastInit, true);
		ref.FAST_Test(im, fast_reference, fastInit, myVector);
		int tempx = -8;
		int curCount = 0;
		for (int i = 0; i < fast_reference.size(); i++)
		{
			KeyPoint kp_opencv = fast_reference[i];
			// 对 kp_opencv 进行操作
			if ((static_cast<int>(kp_opencv.pt.x) / 8) - tempx < time && (static_cast<int>(kp_opencv.pt.x) / 8) >= tempx) {//如果设置8个周期获得一次特征点那么。
				curCount += time - ((static_cast<int>(kp_opencv.pt.x) / 8) - tempx);
				//如果两个特征点的距离小于8那么无事发生，如果大于8那么count需要增加
			}
			tempx = (static_cast<int>(kp_opencv.pt.x) / 8);
		}
		printf("curCount is %d\n", curCount);
		//Mat disimg_opencv;
		//drawKeypoints(im, fast_reference, disimg_opencv);

		count += curCount;

		cv::imshow("Image", im);

		// 等待一段时间
		cv::waitKey(1);
	}
	printf("Count mean is %d\n", count / vstrImageFilenames.size());
	for (int i = 0; i < 9; i++) {
		cout << "8 col sum:" << i << "is count:" << myVector[i] << " data" << myVector[i] / vstrImageFilenames.size() << endl;
	}
	int temp = 0;
	for (int i = 0; i < 9; i++) {
		temp += (myVector[i] / vstrImageFilenames.size()) * i * 8;
	}
	cout << "这是我希望比较小的" <<temp << endl;
}
//不带nms的特征检测
void Verification::FAST_Detection_Test(InputArray _img, int threshold) {
	vector<int> input;
	vector<int> output;
	Verification ver;
	Reference ref;
	ref.FAST_Detection(_img, input, threshold, output);
	String FastPath = "FastDetection\\";
	ver.Output_Data_File(input, FastPath + "ReferenceDataIn.txt", 17);
	ver.Output_Data_File(output, FastPath + "ReferenceDataOut.txt", 1);

}
//得分点的检测
void Verification::CornerScore_Test(InputArray _img, int threshold) {
	vector<int> input;
	vector<int> output;
	Verification ver;
	Reference ref;
	ref.CornerScore(_img, input, threshold, output);
	String FastPath = "CornerScore\\";
	ver.Output_Data_File(input, FastPath + "ReferenceDataIn.txt", 18);
	ver.Output_Data_File(output, FastPath + "ReferenceDataOut.txt", 1);
}
//nms特征筛选的检测
void Verification::NMS_Test(InputArray _img, int threshold) {
	vector<int> input;
	vector<int> output;
	Verification ver;
	Reference ref;
	ref.NMS(_img, input, threshold, output);
	String FastPath = "NMS\\";
	ver.Output_Data_File(input, FastPath + "ReferenceDataIn.txt", 8);
	ver.Output_Data_File(output, FastPath + "ReferenceDataOut.txt", 5);
}

void Verification::FAST_Test(Mat _img, int threshold) {
	vector<int> input;
	vector<int> output;
	vector<int> input2;
	Verification ver;
	Reference ref;
	ref.NMS(_img, input, threshold, output);//FAST得到的结果最后需要经过NMS，使用NMS之后的结果直接输出。
	String FastPath = "Fast\\";
	//ver.Output_Data_File(input, FastPath + "ReferenceDataIn.txt", 8);
	//已经注释了
	//ver.Output_Data_File(output, FastPath + "ReferenceDataOut.txt", 5);
	//ref.ReflectionFillWindow(_img, input2);
	//ver.Output_Data_File(input2, FastPath + "ReferenceDataIn.txt", 8);
}

void Verification::DataGenerateIC_Angle(Mat& image, Point2f pt, std::vector<int>& DataInput) {
	int HALF_PATCH_SIZE = 15;
	vector<int> umax;
	//首先求出圆域
	umax.resize(HALF_PATCH_SIZE + 1);

	int v, v0, vmax = cvFloor(HALF_PATCH_SIZE * sqrt(2.f) / 2 + 1);
	int vmin = cvCeil(HALF_PATCH_SIZE * sqrt(2.f) / 2);
	const double hp2 = HALF_PATCH_SIZE * HALF_PATCH_SIZE;
	for (v = 0; v <= vmax; ++v)
		umax[v] = cvRound(sqrt(hp2 - v * v));

	for (v = HALF_PATCH_SIZE, v0 = 0; v >= vmin; --v)
	{
		while (umax[v0] == umax[v0 + 1])
			++v0;
		umax[v] = v0;
		++v0;
	}
	//从上到下的圆域，pt是中点，通过umax的值选定所需要的数据
	if (pt.x < 15 || pt.y <15 || pt.x > image.rows - 15 || pt.x > image.cols - 15) {
		printf("DataGenerateIC_Angle px error");
		return;
	}
	const uchar* center = &image.at<uchar>(cvRound(pt.y), cvRound(pt.x));
	int step = (int)image.step1();
	for (int y = -HALF_PATCH_SIZE; y <= HALF_PATCH_SIZE; y++) {//根据这个值判断
		for (int x = -HALF_PATCH_SIZE; x <= HALF_PATCH_SIZE; x++) {//如果不在这个范围内就给
			int absy = abs(y);
			int absx = abs(x);
			if (umax[absy] < absx) {
				DataInput.push_back(0);
			}
			else{
				DataInput.push_back(center[x + y * step]);
			}
			printf("y:%d,x%d,%d\n",y,x,DataInput.back());
		}
	}
	for(int i = 0;i < HALF_PATCH_SIZE + 1;i++){
		printf("umax i:%d, umax:%d\n",i ,umax[i]);
	}
}
//随机生成n个数，用于计算灰度质心的mx和my
void Verification::IC_Angle_Test(Mat img, int n) {
	vector<int> input;
	vector<int> output;
	Verification ver;
	Reference ref;
	srand((unsigned)time(NULL));
	cout << "rand num is" << n << endl;
	for (int i = 0; i < n; i++) {
		int mx, my;
		Point2f pt;
		pt.x = (rand() % (img.cols - 40)) + 20;//随机选择n个点
		pt.y = (rand() % (img.rows - 40)) + 20;
		cout << "rand x is" << pt.x << " rand y is" << pt.y << endl;
		//保存输入
		ver.DataGenerateIC_Angle(img, pt, input);

		ref.IC_Angle(img, pt, &mx, &my);//FAST得到的结果最后需要经过NMS，使用NMS之后的结果直接输出。
		output.push_back(my);
		output.push_back(mx);//得到的结果存回
	}
	
	String FastPath = "ICAngle\\";
	//ver.Output_Data_File(input, FastPath + "ReferenceDataIn.txt", 8);
	ver.Output_Data_File(output, FastPath + "ReferenceDataOut.txt", 2, 8);
	ver.Output_Data_File(input, FastPath + "ReferenceDataIn.txt", 31);
}


void Verification::BRIEF_Test(Mat& img, int n) {
	vector<int> input;
	vector<int> output;
	vector<int> res;
	Verification ver;
	Reference ref;
	srand((unsigned)time(NULL));
	cout << "rand num is" << n << endl;
	for (int i = 0; i < n; i++) {
		int mx, my;
		Point2f pt;
		pt.x = (rand() % (img.cols - 40)) + 20;//随机选择n个点
		pt.y = (rand() % (img.rows - 40)) + 20;
		cout << "rand x is" << pt.x << " rand y is" << pt.y << endl;
		//保存输入
		ver.DataGenerateIC_Angle(img, pt, input);

		ref.BRIEF(img, pt, output);//FAST得到的结果最后需要经过NMS，使用NMS之后的结果直接输出。
	}

	String FastPath = "BRIEF\\";
	//ver.Output_Data_File(input, FastPath + "ReferenceDataIn.txt", 8);
	ver.Output_Data_File(output, FastPath + "ReferenceDataOut.txt", 4);
	ver.Output_Data_File(input, FastPath + "ReferenceDataIn.txt", 31);
}



void Verification::tan2_Test(int n) {
	vector<int> input;
	vector<int> output;
	vector<int> res;
	Verification ver;
	Reference ref;
	srand((unsigned)time(NULL));
	cout << "rand num is" << n << endl;
	int count = 0;
	for (int i = 0; i < n; i++) {
		int myres,opencvres, myres1;
		// int x = (rand() % 1024) - 512;//随机选择n个点
		// int y = (rand() % 1024) - 512;
		int x1 = (rand() % (1 << 20) - (1<<19));
		int y1 = (rand() % (1 << 20) - (1<<19));
		//cout << hex << x1 << endl;
		//cout << hex << y1 << endl;
		input.push_back((y1 & 0xFFFFF));
		input.push_back((x1 & 0xFFFFF));
		printf("tan2:%d,%x:,%d,%x:\n",y1,(y1 & 0xFFFFF),x1,(x1 & 0xFFFFF));
		//ref.tan2(y, x, myres);
		ref.tan2_q(y1, x1, myres, 14);
		ref.tan2(y1, x1, myres1);
		float angle = cv::fastAtan2(y1, x1);
		opencvres = static_cast<int>(floorf((angle / 11.25)));
		output.push_back(myres);
		if (myres != opencvres) {
			count++;
			cout << "x:" << x1 << " y:" << y1 << endl; 
			printf("my res:%d,opencv res:%d\n", myres, opencvres);
			cout << myres1 << endl;
			cout << "true valud is" << angle <<"  error count is"<<count<< endl;
			//break;
		}
	}

	String FastPath = "tan2/";
	// //ver.Output_Data_File(input, FastPath + "ReferenceDataIn.txt", 8);
	ver.Output_Data_File(output, FastPath + "ReferenceDataOut.txt", 1);
	ver.Output_Data_File(input, FastPath + "ReferenceDataIn.txt", 2, 5);
}

void Verification::rotate_Test(int n) {
	Reference ref;
	float angle = n * 11.25;
	for (int i = 0; i < 256*4; i+=2) {
		ref.rotate(angle, i);
	}
}

void Verification::rotate_Test1(int n) {//随机生成100个数据 * 32每种情况全部测试一遍。
	Reference ref;
	Verification ver;
	srand((unsigned)time(NULL));
	vector<int> input,output;
	for (int i = 0; i < 32; i++) {
		for (int j = 0; j < n; j++) {
			int a[32];
			for (int k = 0; k < 32; k++) {
				a[k] = (rand() % 128);
				input.push_back(a[k]);
			}
			for (int k = 0; k < 32; k++) {
				output.push_back(a[(k+i)%32]);
			}
		}
	}
	String FastPath = "rotate\\";
	ver.Output_Data_File(output, FastPath + "ReferenceDataOut.txt", 8);
	ver.Output_Data_File(input, FastPath + "ReferenceDataIn.txt", 4);
}



void Verification::RSBrief_Test(Mat& image, std::vector<KeyPoint>& keypoints, vector<vector<uint32_t>> &descriptors) {//给定两个输入分别
	String FastPath = "RSBrief//";
	//产生图片输入
	if (image.cols % 8 != 0) {//必须是8的倍数
		cout << "ReflectionFillWindow error" << endl;
		return;
	}
	Reference ref;
	std::vector<KeyPoint> newKeypoints;

	const int half_patch_size = 8;
	const int half_boundary = 16;
	int bad_points = 0;
	for (auto &kp: keypoints) {
		if (kp.pt.x < half_boundary || kp.pt.y < half_boundary ||
			kp.pt.x >= image.cols - half_boundary || kp.pt.y >= image.rows - half_boundary) {
				continue;
		}
		newKeypoints.push_back(kp);
	}

	ref.RSBrief(image, newKeypoints, descriptors);
	//1、输出图片
	std::vector<int> inputImage = matToVec(image);
	Output_Data_File(inputImage, FastPath + "ReferenceDataInImage.coe", 8);
	//2、输出特征点
	output_FeaturePoint_File(newKeypoints, FastPath + "ReferenceDataInKeypoints.coe");
	//3、输出结果
	std::vector<int> outputDescr;
	for (auto &kp: descriptors) {
		for (auto &temp: kp) {
			outputDescr.push_back(temp);
		}
	}
	Output_Data_File(outputDescr, FastPath + "ReferenceDataOutputDescriptors.coe", 2, 8);
}

void Verification::Top_orb(Mat& image) {//给定两个输入分别
    String FastPath = "Top_orb//";
    //产生图片输入
    Reference ref;
    std::vector<KeyPoint> newKeypoints;

    std::vector<KeyPoint> fast_reference;
    int fastInit = 20;
    std::vector<int> myVector(9);
    ref.FAST_Test(image, fast_reference, fastInit, myVector);
    cout << "fast:opencv" << endl;

    Mat disimg_opencv;
    drawKeypoints(image, fast_reference, disimg_opencv);
    imshow("a", disimg_opencv);
    //fast特征点检测
    const int half_patch_size = 8;
    const int half_boundary = 16;
    int bad_points = 0;
    for (auto &kp: fast_reference) {
        if (kp.pt.x < half_boundary || kp.pt.y < half_boundary ||
            kp.pt.x >= image.cols - half_boundary || kp.pt.y >= image.rows - half_boundary) {
            continue;
        }
        newKeypoints.push_back(kp);
    }

	ofstream ofs1(DIRPATH + FastPath + "OpencvDataOut.txt");
	for (auto kp : newKeypoints) {
		ofs1 << "x:" << kp.pt.x << " y:" << kp.pt.y << " sorce:" << kp.response << endl;
	}

    //删除边缘数据
	Mat cv_res;
	Mat my_res;
	ref.GaussianBlur(image, cv_res, my_res);

    std::vector<std::vector<uint32_t>> descriptors;
    ref.RSBrief(my_res, newKeypoints, descriptors);
    //1、输出图片
    std::vector<int> inputImage = matToVecFill(image);
    Output_Data_File(inputImage, FastPath + "ReferenceDataInImage.coe", 8);
    //2、输出特征点
    output_FeaturePoint_File(newKeypoints, FastPath + "ReferenceDataInKeypoints.coe");
    //3、输出结果
    std::vector<int> outputDescr;
    for (auto &kp: descriptors) {
        for (auto &temp: kp) {
            outputDescr.push_back(temp);
        }
    }
    Output_Data_File(outputDescr, FastPath + "ReferenceDataOutputDescriptors.coe", 2, 8);
}
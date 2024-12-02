#include "Reference.h"
#include "Verification.h"



using namespace std;
/*
	传入原始图像数据和目标图像，还有目标图像的大小。
	要求是灰度图
*/

static int pattern[256 * 4] =
{
	-5,4, -5,13,
	0,11, 4,1,
	-6,-2, -4,-1,
	-11,-4, -3,0,
	4,8, -5,2,
	5,3, -1,12,
	0,7, -5,5,
	9,-2, 8,9,
	-4,5, -2,14,
	2,11, 4,0,
	-6,-1, -4,0,
	-12,-2, -3,1,
	5,7, -5,3,
	5,2, 1,12,
	1,7, -4,6,
	8,-4, 10,7,
	-3,6, 0,14,
	4,10, 4,-1,
	-6,0, -4,1,
	-12,1, -3,1,
	7,6, -4,4,
	6,1, 4,11,
	3,6, -3,7,
	8,-5, 11,5,
	-2,6, 3,14,
	6,9, 4,-1,
	-6,2, -4,1,
	-11,3, -2,2,
	8,4, -3,4,
	6,0, 6,11,
	4,6, -1,7,
	6,-7, 12,3,
	-1,6, 6,13,
	8,8, 4,-2,
	-6,3, -4,2,
	-11,5, -2,2,
	8,3, -2,5,
	6,-1, 8,9,
	5,5, 0,7,
	5,-8, 12,1,
	1,6, 8,11,
	9,6, 3,-3,
	-5,4, -3,3,
	-9,7, -2,2,
	9,1, -1,5,
	5,-2, 9,7,
	6,4, 1,7,
	3,-9, 12,-2,
	2,6, 10,10,
	10,4, 2,-3,
	-4,5, -2,3,
	-8,9, -1,3,
	9,-1, 0,5,
	5,-3, 11,6,
	6,3, 3,7,
	2,-9, 11,-4,
	3,6, 12,7,
	11,2, 2,-4,
	-3,5, -2,4,
	-6,10, -1,3,
	9,-2, 1,5,
	4,-4, 12,3,
	7,1, 4,6,
	0,-9, 10,-6,
	4,5, 13,5,
	11,0, 1,-4,
	-2,6, -1,4,
	-4,11, 0,3,
	8,-4, 2,5,
	3,-5, 12,1,
	7,0, 5,5,
	-2,-9, 9,-8,
	5,4, 14,2,
	11,-2, 0,-4,
	-1,6, 0,4,
	-2,12, 1,3,
	7,-5, 3,5,
	2,-5, 12,-1,
	7,-1, 6,4,
	-4,-8, 7,-10,
	6,3, 14,0,
	10,-4, -1,-4,
	0,6, 1,4,
	1,12, 1,3,
	6,-7, 4,4,
	1,-6, 11,-4,
	6,-3, 7,3,
	-5,-8, 5,-11,
	6,2, 14,-3,
	9,-6, -1,-4,
	2,6, 1,4,
	3,11, 2,2,
	4,-8, 4,3,
	0,-6, 11,-6,
	6,-4, 7,1,
	-7,-6, 3,-12,
	6,1, 13,-6,
	8,-8, -2,-4,
	3,6, 2,4,
	5,11, 2,2,
	3,-8, 5,2,
	-1,-6, 9,-8,
	5,-5, 7,0,
	-8,-5, 1,-12,
	6,-1, 11,-8,
	6,-9, -3,-3,
	4,5, 3,3,
	7,9, 2,2,
	1,-9, 5,1,
	-2,-5, 7,-9,
	4,-6, 7,-1,
	-9,-3, -2,-12,
	6,-2, 10,-10,
	4,-10, -3,-2,
	5,4, 3,2,
	9,8, 3,1,
	-1,-9, 5,0,
	-3,-5, 6,-11,
	3,-6, 7,-3,
	-9,-2, -4,-11,
	6,-3, 7,-12,
	2,-11, -4,-2,
	5,3, 4,2,
	10,6, 3,1,
	-2,-9, 5,-1,
	-4,-4, 3,-12,
	1,-7, 6,-4,
	-9,0, -6,-10,
	5,-4, 5,-13,
	0,-11, -4,-1,
	6,2, 4,1,
	11,4, 3,0,
	-4,-8, 5,-2,
	-5,-3, 1,-12,
	0,-7, 5,-5,
	-9,2, -8,-9,
	4,-5, 2,-14,
	-2,-11, -4,0,
	6,1, 4,0,
	12,2, 3,-1,
	-5,-7, 5,-3,
	-5,-2, -1,-12,
	-1,-7, 4,-6,
	-8,4, -10,-7,
	3,-6, 0,-14,
	-4,-10, -4,1,
	6,0, 4,-1,
	12,-1, 3,-1,
	-7,-6, 4,-4,
	-6,-1, -4,-11,
	-3,-6, 3,-7,
	-8,5, -11,-5,
	2,-6, -3,-14,
	-6,-9, -4,1,
	6,-2, 4,-1,
	11,-3, 2,-2,
	-8,-4, 3,-4,
	-6,0, -6,-11,
	-4,-6, 1,-7,
	-6,7, -12,-3,
	1,-6, -6,-13,
	-8,-8, -4,2,
	6,-3, 4,-2,
	11,-5, 2,-2,
	-8,-3, 2,-5,
	-6,1, -8,-9,
	-5,-5, 0,-7,
	-5,8, -12,-1,
	-1,-6, -8,-11,
	-9,-6, -3,3,
	5,-4, 3,-3,
	9,-7, 2,-2,
	-9,-1, 1,-5,
	-5,2, -9,-7,
	-6,-4, -1,-7,
	-3,9, -12,2,
	-2,-6, -10,-10,
	-10,-4, -2,3,
	4,-5, 2,-3,
	8,-9, 1,-3,
	-9,1, 0,-5,
	-5,3, -11,-6,
	-6,-3, -3,-7,
	-2,9, -11,4,
	-3,-6, -12,-7,
	-11,-2, -2,4,
	3,-5, 2,-4,
	6,-10, 1,-3,
	-9,2, -1,-5,
	-4,4, -12,-3,
	-7,-1, -4,-6,
	0,9, -10,6,
	-4,-5, -13,-5,
	-11,0, -1,4,
	2,-6, 1,-4,
	4,-11, 0,-3,
	-8,4, -2,-5,
	-3,5, -12,-1,
	-7,0, -5,-5,
	2,9, -9,8,
	-5,-4, -14,-2,
	-11,2, 0,4,
	1,-6, 0,-4,
	2,-12, -1,-3,
	-7,5, -3,-5,
	-2,5, -12,1,
	-7,1, -6,-4,
	4,8, -7,10,
	-6,-3, -14,0,
	-10,4, 1,4,
	0,-6, -1,-4,
	-1,-12, -1,-3,
	-6,7, -4,-4,
	-1,6, -11,4,
	-6,3, -7,-3,
	5,8, -5,11,
	-6,-2, -14,3,
	-9,6, 1,4,
	-2,-6, -1,-4,
	-3,-11, -2,-2,
	-4,8, -4,-3,
	0,6, -11,6,
	-6,4, -7,-1,
	7,6, -3,12,
	-6,-1, -13,6,
	-8,8, 2,4,
	-3,-6, -2,-4,
	-5,-11, -2,-2,
	-3,8, -5,-2,
	1,6, -9,8,
	-5,5, -7,0,
	8,5, -1,12,
	-6,1, -11,8,
	-6,9, 3,3,
	-4,-5, -3,-3,
	-7,-9, -2,-2,
	-1,9, -5,-1,
	2,5, -7,9,
	-4,6, -7,1,
	9,3, 2,12,
	-6,2, -10,10,
	-4,10, 3,2,
	-5,-4, -3,-2,
	-9,-8, -3,-1,
	1,9, -5,0,
	3,5, -6,11,
	-3,6, -7,3,
	9,2, 4,11,
	-6,3, -7,12,
	-2,11, 4,2,
	-5,-3, -4,-2,
	-10,-6, -3,-1,
	2,9, -5,1,
	4,4, -3,12,
	-1,7, -6,4,
	9,0, 6,10
};


const float factorPI = (float)(CV_PI / 180.f);
//转为弧度制
void Reference::rotate(float angle, int idx) {//对pattern中的第idx个坐标进行旋转
	float angle1 = (float)angle*factorPI;
	float a = (float)cos(angle1), b = (float)sin(angle1);

	int newy = cvRound(pattern[idx] * b + pattern[idx + 1] * a);
	int newx = cvRound(pattern[idx] * a - pattern[idx + 1] * b);
	cout << idx << " " << newx << " " << newy << endl;
}

//双线性差值的原始代码
void Reference::Resize_Linear(Mat &matSrc, Mat &matDst) {
	//获得目标图像数据
	uchar* dataDst = matDst.data;
	uchar* dataSrc = matSrc.data;
	//获得目标图像一行的字节数
	size_t stepDst = matDst.step;
	size_t stepSrc = matSrc.step;
	//获得原始图像的col和row
	int WidthSrc = matSrc.cols;
	int HiehgtSrc = matSrc.rows;
	//计算缩放比例
	double scale_y = static_cast<double>(matSrc.rows) / matDst.rows;
	double scale_x = static_cast<double>(matSrc.cols) / matDst.cols;
    cout <<"scale_y" << scale_y << endl;
    cout <<"scale_x" << scale_x << endl;
	cout << scale_y * 2048 << endl;
	cout << scale_x * 2048 << endl;
	//开始循环计算
	int a0 = 0;
	int a1 = 0;
	int b0 = 0;
	int b1 = 0;
	for (int j = 0; j < matDst.rows; ++j){
		//i，j是目标图像的像素坐标
		//fy和fx求得目标图像对应原始图像的位置
		float fy = (float)((j + 0.5) * scale_y - 0.5);
		//求完目标距离之后求出距离最近的四个点的相对距离。
		int sy = cvFloor(fy);
		fy -= sy;
		//求得结果应该在[0,1)之间
		//1.2倍的话完全用不到
		if (sy < 0) { //如果映射的目标位置，小于最上侧那么映射的结果只与最上测相关。
			fy = 0, sy = 0;
		}
		if (sy >= HiehgtSrc - 1) {//如果映射的目标位置，大于最下侧那么映射的结果只与最下测相关。
			fy = 1, sy = HiehgtSrc - 2;
		}

		short cbufy[2];
		cbufy[0] = cv::saturate_cast<short>((1.f - fy) * 2048);
		cbufy[1] = 2048 - cbufy[0];
		if (j < 100) {
			cout << "i:" << j << " cbufy0:" << cbufy[0] << endl;
			cout << "i:" << j << " cbufy1:" << cbufy[1] << endl;
			cout << cbufy[0] - a0 << endl;
			cout << cbufy[1] - a1 << endl;
		}
		a0 = cbufy[0];
		a1 = cbufy[1];
		//求的y和（1-y）的数值。
		for (int i = 0; i < matDst.cols; ++i){
			//i，j是目标图像的像素坐标
			//fy和fx求得目标图像对应原始图像的位置
			float fx = (float)((i + 0.5) * scale_x - 0.5);
			int sx = cvFloor(fx);
			fx -= sx;
			//求得结果应该在[0,1)之间
			//对结果的边界进行检查
			if (sx < 0) {
				fx = 0, sx = 0;
			}
			if (sx >= WidthSrc - 1) {
				fx = 1, sx = WidthSrc - 2;
			}

			short cbufx[2];
			cbufx[0] = cv::saturate_cast<short>((1.f - fx) * 2048);
			cbufx[1] = 2048 - cbufx[0];

			if (i < 100 && j == 100) {
				cout << "i:" << j << " cbufx0:" << cbufx[0] << endl;
				cout << "i:" << j << " cbufx1:" << cbufx[1] << endl;
				cout << cbufx[0] - b0 << endl;
				cout << cbufx[1] - b1 << endl;
			}
			b0 = cbufx[0];
			b1 = cbufx[1];

			//求得目标像素跟原始像素的相对距离
			//计算结果
			for (int k = 0; k < matSrc.channels(); ++k) {
				//				*(dataDst + j * stepDst + matSrc.channels() * i + k) = round((double)(
				//					*(dataSrc + sy * stepSrc + matSrc.channels() * sx + k) * cbufx[0] * cbufy[0] +
				//					*(dataSrc + (sy + 1)*stepSrc + matSrc.channels() * sx + k) * cbufx[0] * cbufy[1] +
				//					*(dataSrc + sy * stepSrc + matSrc.channels() * (sx + 1) + k) * cbufx[1] * cbufy[0] +
				//					*(dataSrc + (sy + 1)*stepSrc + matSrc.channels() * (sx + 1) + k) * cbufx[1] * cbufy[1]) / (1 << 22));
				//				*(dataDst + j * stepDst + matSrc.channels() * i + k) = (
				//					*(dataSrc + sy * stepSrc + matSrc.channels() * sx + k) * cbufx[0] * cbufy[0] +
				//					*(dataSrc + (sy + 1)*stepSrc + matSrc.channels() * sx + k) * cbufx[0] * cbufy[1] +
				//					*(dataSrc + sy * stepSrc + matSrc.channels() * (sx + 1) + k) * cbufx[1] * cbufy[0] +
				//					*(dataSrc + (sy + 1)*stepSrc + matSrc.channels() * (sx + 1) + k) * cbufx[1] * cbufy[1]) >> 22;
				*(dataDst + j * stepDst + matSrc.channels() * i + k) = (
					(((((*(dataSrc + sy * stepSrc + matSrc.channels() * sx + k) * cbufx[0] +
						*(dataSrc + sy * stepSrc + matSrc.channels() * (sx + 1) + k) * cbufx[1]) >> 4) * cbufy[0]) >> 16) +
						((((*(dataSrc + (sy + 1)*stepSrc + matSrc.channels() * sx + k) * cbufx[0] +
							*(dataSrc + (sy + 1)*stepSrc + matSrc.channels() * (sx + 1) + k) * cbufx[1]) >> 4) * cbufy[1]) >> 16) + 2) >> 2
					);
			}
		}
	}
}

//测试1：将Resize_Linear中涉及的所有数全部定点化，为了简单规定所采用的数需要有负数，因此采用32位的有符号数
//最高位为符号位，scale_y = 这个定点化，包含了多个多个缩放纬度的参数，比较通用。
void Reference::Resize_Linear_Test(Mat &matSrc, Mat &matDst) {
	//获得目标图像数据
	int offset = 11;
	uchar* dataDst = matDst.data;
	uchar* dataSrc = matSrc.data;
	//获得目标图像一行的字节数
	int stepDst = matDst.step;
	int stepSrc = matSrc.step;
	//获得原始图像的col和row
	int WidthSrc = matSrc.cols;
	int HiehgtSrc = matSrc.rows;
	//计算缩放比例
	//修改1，在这里进行修改，缩放的话这里的位宽一定会缩小，因此采用26位之内的位宽即可，这里可以使用26 + 11位的位宽来使用。37 + 26来计算。最高位采用符号位因此使用38 + 26来计算
	long long scale_y = (long long)(((double)((long long)matSrc.rows << offset) / matDst.rows) + 0.5);
	long long scale_x = (long long)(((double)((long long)matSrc.cols << offset) / matDst.cols) + 0.5);//左移20位
	//int scale_y = cvRound(double(matSrc.rows << 20) / matDst.rows);
	//int scale_x = cvRound(double(matSrc.cols << 20) / matDst.cols);//左移20位
	cout <<"scale_y" << scale_y << endl;
	cout <<"scale_x" << scale_x << endl;
	//开始循环计算
	for (int j = 0; j < matDst.rows; ++j) {
		//i，j是目标图像的像素坐标
		//fy和fx求得目标图像对应原始图像的位置
		//修改2，对0.5进行修改，左移1位
		//1-----------------------求出fy
		long long fy = ((2 * j + 1) * scale_y - (1l << offset)); //一共左移了21位
		//求完目标距离之后求出距离最近的四个点的相对距离。
		//修改后第24位是小数部分，因此直接取值即可。
		int sy = fy >> (offset + 1);//取整数部分。
		fy = (fy % (1 << (offset + 1)));//取小数部分。
		//求得结果应该在[0,1)之间
		//2-----------------------求出fy，根据这里的值判断是否跳过，如果不相等则跳过这里。
		if (sy < 0) { //如果映射的目标位置，小于最上侧那么映射的结果只与最上测相关。
			fy = 0, sy = 0;
		}
		if (sy >= HiehgtSrc - 1) {//如果映射的目标位置，大于最下侧那么映射的结果只与最下测相关。
			fy = 1 << (offset + 1) , sy = HiehgtSrc - 2; //这里时固定值
		}
		short cbufy[2];
		//cbufy[0] = cv::saturate_cast<short>((2048 - round((double)fy / (1 << (offset - 10)))));//这里修改为+1取整
		//cbufy[1] = 2048 - cbufy[0];
		cbufy[0] = (short)((4096 - (fy >> (offset - 11))) >> 1);
		cbufy[1] = ((fy >> (offset - 11)) + 1) >> 1;
		if (j < 100) {
			cout << "i:" << j << " cbufx0:" << cbufy[0] << endl;
			cout << "i:" << j << " cbufx1:" << cbufy[1] << endl;
		}
		//3---------------------将上一个数据寄存一个周期，然后开始计算，应该
		//求的y和（1-y）的数值。
		for (int i = 0; i < matDst.cols; ++i) {
			//i，j是目标图像的像素坐标
			//fy和fx求得目标图像对应原始图像的位置
			long long fx = ((2 * i + 1) * scale_x - (1 << offset));
			int sx = fx >> (offset + 1);//取整数部分
			fx = (fx % (1 << (offset + 1)));//取小数部分
			//求得结果应该在[0,1)之间
			//对结果的边界进行检查
			if (sx < 0) {
				fx = 0, sx = 0;
			}
			if (sx >= WidthSrc - 1) {
				fx = 1 << (offset + 1), sx = WidthSrc - 2;
			}
			//求得结果之后直接给出，2级流水给出，因为不做图像放大，所以这个数据一定不会再次利用到，直接丢弃。而且需要做一个判断判断接收到的数据是否在当前计算的范围内。
			//如果不在需要，无视数据。
			//要求确保下一级模块的ready一直拉高，上级模块的valid值允许拉低，一旦下级模块ready信号拉低会造成数据丢失，所以采用fifo的形式，
			//采用例如64个数据存储的形式可以一次性保存较多的数据。
			short cbufx[2];
			//cbufx[0] = cv::saturate_cast<short>((2048 - round((double)fx / (1<<(offset - 10)))));
			//cbufx[1] = 2048 - cbufx[0];
			//cbufx[0] = cv::saturate_cast<short>((4096 - (fx >> (offset - 11))) >> 1);
			//cbufx[1] = ((fx >> (offset - 11)) + 1) >> 1;
			cbufx[0] = (short)(2048 - (fx >> (offset - 10)));
			cbufx[1] = (fx >> (offset - 10));
			//求得目标像素跟原始像素的相对距离
			//计算结果
			for (int k = 0; k < matSrc.channels(); ++k) {
				//				*(dataDst + j * stepDst + matSrc.channels() * i + k) = round((double)(
				//					*(dataSrc + sy * stepSrc + matSrc.channels() * sx + k) * cbufx[0] * cbufy[0] +
				//					*(dataSrc + (sy + 1)*stepSrc + matSrc.channels() * sx + k) * cbufx[0] * cbufy[1] +
				//					*(dataSrc + sy * stepSrc + matSrc.channels() * (sx + 1) + k) * cbufx[1] * cbufy[0] +
				//					*(dataSrc + (sy + 1)*stepSrc + matSrc.channels() * (sx + 1) + k) * cbufx[1] * cbufy[1]) / (1 << 22));
				//				*(dataDst + j * stepDst + matSrc.channels() * i + k) = (
				//					*(dataSrc + sy * stepSrc + matSrc.channels() * sx + k) * cbufx[0] * cbufy[0] +
				//					*(dataSrc + (sy + 1)*stepSrc + matSrc.channels() * sx + k) * cbufx[0] * cbufy[1] +
				//					*(dataSrc + sy * stepSrc + matSrc.channels() * (sx + 1) + k) * cbufx[1] * cbufy[0] +
				//					*(dataSrc + (sy + 1)*stepSrc + matSrc.channels() * (sx + 1) + k) * cbufx[1] * cbufy[1]) >> 22;
				*(dataDst + j * stepDst + matSrc.channels() * i + k) = (
					(((((*(dataSrc + sy * stepSrc + matSrc.channels() * sx + k) * cbufx[0] +
						*(dataSrc + sy * stepSrc + matSrc.channels() * (sx + 1) + k) * cbufx[1]) >> 4) * cbufy[0]) >> 16) +
						((((*(dataSrc + (sy + 1)*stepSrc + matSrc.channels() * sx + k) * cbufx[0] +
							*(dataSrc + (sy + 1)*stepSrc + matSrc.channels() * (sx + 1) + k) * cbufx[1]) >> 4) * cbufy[1]) >> 16) + 2) >> 2
					);
			}
		}
	}
}
//特意优化的1.25倍缩放，量化参数少
void Reference::Resize_Linear_Test_1_2(Mat &matSrc, Mat &matDst) {
	//获得目标图像数据
	int offset = 2;
	uchar* dataDst = matDst.data;
	uchar* dataSrc = matSrc.data;
	//获得目标图像一行的字节数
	int stepDst = matDst.step;
	int stepSrc = matSrc.step;
	//获得原始图像的col和row
	int WidthSrc = matSrc.cols;
	int HiehgtSrc = matSrc.rows;
	//计算缩放比例
	//修改1，在这里进行修改，缩放的话这里的位宽一定会缩小，因此采用26位之内的位宽即可，这里可以使用26 + 11位的位宽来使用。37 + 26来计算。最高位采用符号位因此使用38 + 26来计算
//	long long scale_y = (long long)(((double)((long long)matSrc.rows << offset) / matDst.rows));
//	long long scale_x = (long long)(((double)((long long)matSrc.cols << offset) / matDst.cols));//左移20位
    long long scale_y = 5;
    long long scale_x = 5;
	//int scale_y = cvRound(double(matSrc.rows << 20) / matDst.rows);
	//int scale_x = cvRound(double(matSrc.cols << 20) / matDst.cols);//左移20位
	cout << "scale_y" << scale_y << endl;
	cout << "scale_x" << scale_x << endl;
	//开始循环计算
	for (int j = 0; j < matDst.rows; ++j) {
		//i，j是目标图像的像素坐标
		//fy和fx求得目标图像对应原始图像的位置
		//修改2，对0.5进行修改，左移1位
		//1-----------------------求出fy
		long long fy = ((2 * j + 1) * scale_y - (1l << offset)); //一共左移了21位
		//求完目标距离之后求出距离最近的四个点的相对距离。
		//修改后第24位是小数部分，因此直接取值即可。
		int sy = fy >> (offset + 1);//取整数部分。
		fy = (fy % (1 << (offset + 1)));//取小数部分。
		//求得结果应该在[0,1)之间
		//2-----------------------求出fy，根据这里的值判断是否跳过，如果不相等则跳过这里。
		if (sy < 0) { //如果映射的目标位置，小于最上侧那么映射的结果只与最上测相关。
			fy = 0, sy = 0;
		}
		if (sy >= HiehgtSrc - 1) {//如果映射的目标位置，大于最下侧那么映射的结果只与最下测相关。
			fy = 1 << (offset + 1), sy = HiehgtSrc - 2; //这里时固定值
		}
		short cbufy[2];
		//cbufy[0] = cv::saturate_cast<short>((2048 - round((double)fy / (1 << (offset - 10)))));//这里修改为+1取整
		//cbufy[1] = 2048 - cbufy[0];
		cbufy[0] = (short)((1 << (offset + 1)) - fy) << (10 - offset);
		cbufy[1] = (fy << (10 - offset));
		if (j < 100) {
			cout << "i:" << j << " cbufy0:" << cbufy[0] << endl;
			cout << "i:" << j << " cbufy1:" << cbufy[1] << endl;
		}
		//3---------------------将上一个数据寄存一个周期，然后开始计算，应该
		//求的y和（1-y）的数值。
		for (int i = 0; i < matDst.cols; ++i) {
			//i，j是目标图像的像素坐标
			//fy和fx求得目标图像对应原始图像的位置
			long long fx = ((2 * i + 1) * scale_x - (1l << offset));
			int sx = fx >> (offset + 1);//取整数部分
			fx = (fx % (1 << (offset + 1)));//取小数部分
			//求得结果应该在[0,1)之间
			//对结果的边界进行检查
			if (sx < 0) {
				fx = 0, sx = 0;
			}
			if (sx >= WidthSrc - 1) {
				fx = 1 << (offset + 1), sx = WidthSrc - 2;
			}
			if (i < 100 && j == 100) {
				cout << "i:" << j << " fx:" << fx << endl;
				//cout << "i:" << j << " cbufx1:" << cbufx[1] << endl;
			}
			//求得结果之后直接给出，2级流水给出，因为不做图像放大，所以这个数据一定不会再次利用到，直接丢弃。而且需要做一个判断判断接收到的数据是否在当前计算的范围内。
			//如果不在需要，无视数据。
			//要求确保下一级模块的ready一直拉高，上级模块的valid值允许拉低，一旦下级模块ready信号拉低会造成数据丢失，所以采用fifo的形式，
			//采用例如64个数据存储的形式可以一次性保存较多的数据。
			short cbufx[2];
			//cbufx[0] = cv::saturate_cast<short>((2048 - round((double)fx / (1<<(offset - 10)))));
			//cbufx[1] = 2048 - cbufx[0];
			//cbufx[0] = cv::saturate_cast<short>((4096 - (fx >> (offset - 11))) >> 1);
			//cbufx[1] = ((fx >> (offset - 11)) + 1) >> 1;
			cbufx[0] = (short)(((1 << (offset + 1)) - fx) << (10 - offset));
			cbufx[1] = (fx << (10 - offset));
			//求得目标像素跟原始像素的相对距离
			//计算结果
			for (int k = 0; k < matSrc.channels(); ++k) {
				//				*(dataDst + j * stepDst + matSrc.channels() * i + k) = round((double)(
				//					*(dataSrc + sy * stepSrc + matSrc.channels() * sx + k) * cbufx[0] * cbufy[0] +
				//					*(dataSrc + (sy + 1)*stepSrc + matSrc.channels() * sx + k) * cbufx[0] * cbufy[1] +
				//					*(dataSrc + sy * stepSrc + matSrc.channels() * (sx + 1) + k) * cbufx[1] * cbufy[0] +
				//					*(dataSrc + (sy + 1)*stepSrc + matSrc.channels() * (sx + 1) + k) * cbufx[1] * cbufy[1]) / (1 << 22));
				//				*(dataDst + j * stepDst + matSrc.channels() * i + k) = (
				//					*(dataSrc + sy * stepSrc + matSrc.channels() * sx + k) * cbufx[0] * cbufy[0] +
				//					*(dataSrc + (sy + 1)*stepSrc + matSrc.channels() * sx + k) * cbufx[0] * cbufy[1] +
				//					*(dataSrc + sy * stepSrc + matSrc.channels() * (sx + 1) + k) * cbufx[1] * cbufy[0] +
				//					*(dataSrc + (sy + 1)*stepSrc + matSrc.channels() * (sx + 1) + k) * cbufx[1] * cbufy[1]) >> 22;
				*(dataDst + j * stepDst + matSrc.channels() * i + k) = (//
					(((((*(dataSrc + sy * stepSrc + matSrc.channels() * sx + k) * cbufx[0] +
						*(dataSrc + sy * stepSrc + matSrc.channels() * (sx + 1) + k) * cbufx[1]) >> 4) * cbufy[0]) >> 16) +
						((((*(dataSrc + (sy + 1)*stepSrc + matSrc.channels() * sx + k) * cbufx[0] +
							*(dataSrc + (sy + 1)*stepSrc + matSrc.channels() * (sx + 1) + k) * cbufx[1]) >> 4) * cbufy[1]) >> 16) + 2) >> 2
					);
			}
		}
	}
}

vector<int> Reference::Windows_Test(Mat &matSrc, int size_h, int size_w, int data_num) {
	vector<int> windows;
	if (matSrc.empty() || matSrc.channels() != 1 || matSrc.cols % data_num != 0) {//必须满足一行的个数是data_num的倍数例如8或者10
		cout << "Windows_Test_Error" << endl;
		return windows;
	}
	size_t stepDst = matSrc.step;
	int width = matSrc.cols;
	int height = matSrc.rows;
	uchar* dataSrc = matSrc.data;
	cout << stepDst << endl;
	cout << width << endl;
	for (int h = 0; h <= height - size_h; h++ ) {
		for (int w = 0; w <= width - size_w * data_num; w += data_num) {
			for (int i = 0; i < size_h; i++) {
				for (int j = 0; j < size_w; j++) {
					for (int k = 0; k < data_num; k++) {
						windows.push_back((int)dataSrc[(h + i) * width + w + j * data_num + k]);
					}
				}
			}
		}
	}
	return windows;
}

//vector<int> Reference::Fast_Data_Generater(Mat &matSrc, int data_num) {
//	bool resize = true;
//	vector<int> windows;
//	if (matSrc.empty() || matSrc.channels() != 1 || matSrc.cols % 10 != 0) {
//		cout << "Windows_Test_Error" << endl;
//		return windows;
//	}
//	size_t stepDst = matSrc.step;
//	int width = matSrc.cols;
//	int height = matSrc.rows;
//	uchar* dataSrc = matSrc.data;
//	cout << stepDst << endl;
//	cout << width << endl;
//	for (int h = 0; h <= height - size_h; h++) {
//		if (resize && h % 5 == 4)
//			continue;
//		for (int w = 0; w <= width - size_w * data_num; w += data_num) {
//			for (int i = 0; i < size_h; i++) {
//				for (int j = 0; j < size_w; j++) {
//					for (int k = 0; k < data_num; k++) {
//						windows.push_back((int)dataSrc[(h + i) * width + w + j * data_num + k]);
//					}
//				}
//			}
//		}
//	}
//	return windows;
//}

vector<Mat> Reference::ComputePyramid(cv::Mat image)
{
	float mvInvScaleFactor[2] = { 1 / 1.2f , 1 / 1.44f };
	vector<Mat> mvImagePyramid(2);
	int EDGE_THRESHOLD = 19;
	for (int level = 0; level < 2; ++level)
	{
		float scale = mvInvScaleFactor[level];
		Size sz(cvRound((float)image.cols*scale), cvRound((float)image.rows*scale));
		Size wholeSize(sz.width + EDGE_THRESHOLD * 2, sz.height + EDGE_THRESHOLD * 2);
		Mat temp(wholeSize, image.type()), masktemp;
		mvImagePyramid[level] = temp(Rect(EDGE_THRESHOLD, EDGE_THRESHOLD, sz.width, sz.height));

		// Compute the resized image
		if (level != 0)
		{
			resize(mvImagePyramid[level - 1], mvImagePyramid[level], sz, 0, 0, INTER_LINEAR);

			printf("mvImagePyramid : col:%d,row:%d\n", mvImagePyramid[level].cols, mvImagePyramid[level].rows);

			copyMakeBorder(mvImagePyramid[level], temp, EDGE_THRESHOLD, EDGE_THRESHOLD, EDGE_THRESHOLD, EDGE_THRESHOLD,
				BORDER_REFLECT_101 + BORDER_ISOLATED);

			printf("temp : col:%d,row:%d\n", temp.cols, temp.rows);
			imshow("Test", temp);
		}
		else
		{
			copyMakeBorder(image, temp, EDGE_THRESHOLD, EDGE_THRESHOLD, EDGE_THRESHOLD, EDGE_THRESHOLD,
				BORDER_REFLECT_101);
		}
	}
	return mvImagePyramid;
}


int cornerScore(Mat &img, int i, int j, int threshold) {
	int x[25] = { 0, 1, 2, 3, 3, 3, 2, 1, 0, -1, -2, -3, -3, -3, -2, -1, 0, 1, 2, 3, 3, 3, 2, 1, 0 };
	int y[25] = { 3, 3, 2, 1, 0, -1, -2, -3, -3, -3,-2 ,-1, 0, 1, 2, 3, 3, 3, 2, 1, 0, -1, -2, -3,-3 };
	const int K = 8, N = K * 3 + 1;
	uchar* dataSrc = img.data;
	size_t step = img.step;
	int k, v = *(dataSrc + (i) * step + j);
	short d[N];
	for (k = 0; k < N; k++)    //计算与周围16个像素点的差值,保存在d[k]中  
		d[k] = (short)(v - *(dataSrc + (i + y[k]) * step + j + x[k]));
	//分别计算两
	int a0 = threshold;
	for (k = 0; k < 16; k += 2)  //周围像素小于中心点像素  
	{
		int a = std::min((int)d[k + 1], (int)d[k + 2]);
		a = std::min(a, (int)d[k + 3]);
		if (a <= a0)
			continue;
		a = std::min(a, (int)d[k + 4]);
		a = std::min(a, (int)d[k + 5]);
		a = std::min(a, (int)d[k + 6]);
		a = std::min(a, (int)d[k + 7]);
		a = std::min(a, (int)d[k + 8]);
		a0 = std::max(a0, std::min(a, (int)d[k]));
		a0 = std::max(a0, std::min(a, (int)d[k + 9]));
	}

	int b0 = -a0;
	for (k = 0; k < 16; k += 2)//周围像素点大于中心像素点  
	{
		int b = std::max((int)d[k + 1], (int)d[k + 2]);
		b = std::max(b, (int)d[k + 3]);
		b = std::max(b, (int)d[k + 4]);
		b = std::max(b, (int)d[k + 5]);
		if (b >= b0)
			continue;
		b = std::max(b, (int)d[k + 6]);
		b = std::max(b, (int)d[k + 7]);
		b = std::max(b, (int)d[k + 8]);

		b0 = std::min(b0, std::max(b, (int)d[k]));
		b0 = std::min(b0, std::max(b, (int)d[k + 9]));
	}

	threshold = -b0 - 1;

	return threshold;
}

int GPcountFast = 0;
int GPcountNMSFast = 0;
int GPperform = 0;
int GPMaxCount = 0;
int GPMinCount = 9999999;

int Reference::getPerf(int idx){
    switch (idx)
    {
    	case 0:  return GPcountFast;
    	case 1:  return GPcountNMSFast;
        case 2:   return GPperform;
        case 3:  return GPMaxCount;
        case 4:  return GPMinCount;
        default:   return -1;
    }
}

int Reference::clearPerf(int idx){
    switch (idx)
    {
    	case 0:  GPcountFast = 0; return 1;
    	case 1:  GPcountNMSFast = 0;return 1;
        case 2:  GPperform = 0;return 1;
        case 3:  GPMaxCount = 0;
        case 4:  GPMinCount = 9999999;
        default:   return -1;
    }
}


//首先对fast进行测试，特征检测，筛选，
void Reference::FAST_Test(InputArray _img, std::vector<KeyPoint>& keypoints, int threshold, std::vector<int>& temp)
{

	Mat img = _img.getMat();
	const int K = 8, N = 25;
	int sum[9] = { 0 };
	int groupSum = 0;
	uchar* dataSrc = img.data;
	//获得目标图像一行的字节数
	size_t step = img.step;
	int countFast = 0;
	int countNMSFast = 0;


	int kkk[10] = { 0 };
	threshold = std::min(std::max(threshold, 0), 255);//保证阈值在0-255之间。  
	int x[25] = {0, 1, 2, 3, 3, 3, 2, 1, 0, -1, -2, -3, -3, -3, -2, -1, 0, 1, 2, 3, 3, 3, 2, 1, 0};
	int y[25] = {3, 3, 2, 1, 0, -1, -2, -3, -3, -3,-2 ,-1, 0, 1, 2, 3, 3, 3, 2, 1, 0, -1, -2, -3,-3};
	vector<vector<int>> matrix(img.rows, vector<int>(img.cols, 0));
	for (int i = 3; i < img.rows - 3; i++) {
		for (int j = 3; j < img.cols - 3; j++) {
			int v = *(dataSrc + i * step + j);
			int vt = v - threshold;
			int count = 0;
			for (int k = 0; k < 25; k++) {
				int dark = *(dataSrc + (i + y[k]) * step + j + x[k]);
				if (dark < vt)
				{
					if (++count > K)
					{
						groupSum++;
						matrix[i][j] = cornerScore(img, i, j, threshold);
						break;
					}
				}
				else
					count = 0;
			}


			count = 0;
			vt = v + threshold;
			for (int k = 0; k < 25; k++) {
				int light = *(dataSrc + (i + y[k]) * step + j + x[k]);
				if (light > vt)
				{
					if (++count > K)
					{
						groupSum++;
						matrix[i][j] = cornerScore(img, i, j, threshold);
						break;
					}
				}
				else
					count = 0;
			}
			if (j % 8 == 2 || j == img.cols - 4) {
				sum[groupSum]++;
				groupSum = 0;
			}
		}
	}
	for (int i = 3; i < img.rows - 3; i++) {
		for (int j = 3; j < img.cols - 3; j++) {
			int score = matrix[i][j];
			if (score > 0) {
				countFast++;
			}
			if(score > matrix[i+1][j+1] && score > matrix[i + 1][j] && score > matrix[i + 1][j - 1] && 
				score > matrix[i][j - 1] && score > matrix[i][j + 1] && 
				score > matrix[i - 1][j + 1] && score > matrix[i - 1][j] && score > matrix[i - 1][j - 1] &&
				score > 0){
				keypoints.push_back(KeyPoint((float)j, (float)(i), 7.f, -1, (float)score));
				countNMSFast++;
			}
		}
	}
	cout << "fast num brefor nms:" << countFast << endl;
	cout << "fast num after nms:" << countNMSFast << endl;
    int clk = 0;
	for (int i = 0; i < 9; i++) {
		cout << "8 col sum:" << i << "is count:" << sum[i] << endl;
		temp[i] += sum[i];
		if(i == 0){
		    clk += sum[i];
		}else{
		    clk += sum[i] * i;
            GPperform += sum[i] * (i - 1);
		}
	}
    cout << "fast small time:" << clk << endl;
    cout << "fast full  time:" << img.rows * img.cols / 8 << endl;
    cout << "%" << (double)clk/(img.rows * img.cols / 8)<<endl;

    //-----------------------perf------------------------
    GPcountFast+=countFast;
    GPcountNMSFast+=countNMSFast;
    GPMaxCount = max(GPMaxCount, countNMSFast);
    GPMinCount = min(GPMinCount, countNMSFast);
	return;
}

//padding的对称填充，得到的windows
void Reference::ReflectionFillWindow(Mat &_img, std::vector<int64_t>& keypoints,int h, int w)
{
//	int kernelSize = 7;
	int udRadium = h/2;
	int lfRadium = w/2;
	cv::Mat paddedImg;
	//paddedImg = _img.clone();
	int rows = _img.rows;
	int cols = _img.cols;

	std::vector<std::vector<int64_t>> newVec(rows + udRadium * 2, std::vector<int64_t>(cols/8 + lfRadium * 2));
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j+=8) {
			int64_t a = 0;
			for(int k = 0; k < 8; k++){
				a = (a << 8) + _img.at<uchar>(i, j + 7 - k);
			}
			newVec[i + udRadium][j/8 + lfRadium] = a;
		}
	}

	rows = rows + udRadium * 2;
	cols = cols / 8 + lfRadium * 2;
	//先映射左右
	for (int i = 0; i < rows; i++) {
		for(int j = 0;j < lfRadium;j++){
			newVec[i][j] = newVec[i][j + (lfRadium - j) * 2];
			newVec[i][cols - 1 - j] = newVec[i][cols - 1 - j - (lfRadium - j) * 2];
		}
	}
	//映射上下
	for (int i = 0; i < cols; i++) {
		for(int j = 0;j < udRadium;j++){
			newVec[j][i] = newVec[j + (udRadium - j) * 2][i];
			newVec[rows - j - 1][i] = newVec[rows - j - 1 - (udRadium - j) * 2][i];
		}
	}
	//填充数据
	for (int i = 0; i < rows - udRadium * 2; i++) {
		for (int j = 0; j < cols - lfRadium * 2; j++) {
			for (int a = 0; a < (udRadium *2 + 1); ++a) {
				for (int b = 0; b < (lfRadium * 2 + 1); ++b) {
					keypoints.push_back(newVec[i + a][j + b]);
				}
			}
		}
	}
}

void Reference::noPaddingWindow(Mat &_img, std::vector<int>& keypoints)
{
	int kernelSize = 7;
	int radium = 3;
	int rows = _img.rows;
	int cols = _img.cols;
	for (int i = 0; i < rows - 6; ++i) {
		for (int j = 0; j < cols - 7; j += 8) {
			for (int a = 0; a < 7; ++a) {
				for (int b = 0; b < 8; ++b) {
					keypoints.push_back(_img.at<uchar>(i + a, j + b));
				}
			}
		}
	}
}

void Reference::FastData(Mat &image, std::vector<int>& keypoints)
{
	if (image.cols % 8 != 0) {
		cout << "FastData error" << endl;
	}
	int kernelSize = 7;
	int radium = 3;
	cv::Mat paddedImg;

	cv::copyMakeBorder(image, paddedImg, radium, radium, radium, radium, cv::BORDER_REFLECT_101);
	int rows = paddedImg.rows;
	int cols = paddedImg.cols;
	for (int i = 0; i < rows - 6; ++i) {
		for (int j = 0; j < cols - 7; j += 8) {
			for (int a = 0; a < 7; ++a) {
				for (int b = 0; b < 3 + 8 + 3; ++b) {
					keypoints.push_back(paddedImg.at<uchar>(i + a, j + b));
				}
			}
		}
	}
}


void Reference::GaussianBlur(Mat &image, Mat &my_res) {

//	if (image.cols % 8 != 0) {
//		cout << "GaussianBlur error" << endl;
//	}
    uchar data[] = { 9, 17, 24, 28, 24, 17, 9 };
    Mat mat(7, 1, CV_8U, data);
    Mat outputMat;
    Reference ref;
    Mat paddingMat;
    int kernelSize = 7;
    int radium = 3;
    Mat result2(image.rows, image.cols, CV_8U, Scalar(0));
    cv::copyMakeBorder(image, paddingMat, radium, radium, radium, radium, cv::BORDER_REFLECT_101);

    int rows = paddingMat.rows;
    int cols = paddingMat.cols;

    for (int i = radium; i < rows - radium; ++i) {
        for (int j = radium; j < cols - radium; ++j) {
            int sum[7] = { 0 };
            int countSum = 0;
            for (int a = -radium; a <= radium; ++a) {//这里是列
                for (int b = -radium; b <= radium; ++b) {//这里是行
                    int value1 = paddingMat.at<uchar>(i + b, j + a);
                    int value2 = mat.at<uchar>(b + radium, 0);
                    //cout << a.cols << " " << a.rows << endl;
                    sum[a + radium] += value1 * value2;
                }
                int value1 = sum[a + radium];
                int value2 = mat.at<uchar>(a + radium, 0);
                countSum += value1 * value2;
            }

            double temp = (double)countSum / (128 * 128);
            result2.at<uchar>(i - radium, j - radium) = static_cast<uchar>(cvRound(temp));
        }
    }
    my_res = result2;
}

void Reference::GaussianBlur(Mat &image, Mat &cv_res, Mat &my_res) {
	
//	if (image.cols % 8 != 0) {
//		cout << "GaussianBlur error" << endl;
//	}
	uchar data[] = { 9, 17, 24, 28, 24, 17, 9 };
	Mat mat(7, 1, CV_8U, data);
	Mat outputMat;
	Reference ref;
	Mat paddingMat;
	int kernelSize = 7;
	int radium = 3;
	Mat result2(image.rows, image.cols, CV_8U, Scalar(0));
	cv::copyMakeBorder(image, paddingMat, radium, radium, radium, radium, cv::BORDER_REFLECT_101);
	
	int rows = paddingMat.rows;
	int cols = paddingMat.cols;
	
	for (int i = radium; i < rows - radium; ++i) {
		for (int j = radium; j < cols - radium; ++j) { 
			int sum[7] = { 0 };
			int countSum = 0; 
			for (int a = -radium; a <= radium; ++a) {//这里是列
				for (int b = -radium; b <= radium; ++b) {//这里是行
					int value1 = paddingMat.at<uchar>(i + b, j + a);
					int value2 = mat.at<uchar>(b + radium, 0);
					//cout << a.cols << " " << a.rows << endl;
					sum[a + radium] += value1 * value2;
				}
				int value1 = sum[a + radium];
				int value2 = mat.at<uchar>(a + radium, 0);
				countSum += value1 * value2;
			}

			double temp = (double)countSum / (128 * 128);
			result2.at<uchar>(i - radium, j - radium) = static_cast<uchar>(cvRound(temp));
		}
	}
	imshow("my:GaussianBlur", result2);

	cv::GaussianBlur(image, outputMat, Size(7, 7), 2, 2, BORDER_REFLECT_101);
	cout << outputMat.rows << " " << outputMat.cols << endl;
	imshow("cv:GaussianBlur", outputMat);
	imshow("imgae", image);
	cv_res = outputMat;
	my_res = result2;

	//resize(src, dsrc, Size(640, 640), 0, 0, INTER_LINEAR);

	//	for (int row = 0; row < 7; row++) {
	//		for (int col = 0; col < 7; col++) {
	//			double value0 = dsrc.at<double>(row, 0);
	//			double value1 = dsrc.at<double>(col, 0);
	//			if (row == 0 && col == 0) {
	//				k = value0 * value1;
	//			}
	//			printf("%f ", (value0 * value1) / k);
	//		}
	//		printf("\n");
	//	}
		//vector<Mat> cp = ref.ComputePyramid(dsrc);
		//ver.ComputePyramid(dsrc);
		//for (int i = 0; i < 2; i++) {
		//	printf("col:%d,row:%d\n", cp[i].cols, cp[i].rows);
		//}

		//Mat image_ori = src(rect);
		//Size dsize(640 / 1.25, 640 / 1.25);
		//ver.BitwidthConversion_Test_Verification(dsrc);

		//64 54 45 38 31 26 22 18
		//ver.FAST_Test_Verification(dsrc);
	//	auto arr1 = ver.matToVec(src);
	//	ofstream ofs(DIRPATH + "Resize\\" + "TestDataSrc.txt");
	//	for (int row = 0; row < src.rows; row++) {
	//		for (int col = 0; col < src.cols; col++) {
	//			for (int dim = 0; dim < src.channels(); dim++) {
	//				ofs << hex << setfill('0') << setw(2) << arr1[row * src.cols * src.channels() + col * src.channels() + dim];
	//			}
	//			ofs << endl;
	//		}
	//	}

}


//特征点检测和得分点计算没有nms
void Reference::FAST_Detection(InputArray _img, std::vector<int>& keypoints, int threshold, std::vector<int>& temp)
{

	Mat img = _img.getMat();
	const int K = 8, N = 25;
	int sum[9] = { 0 };
	int groupSum = 0;
	uchar* dataSrc = img.data;
	//获得目标图像一行的字节数
	size_t step = img.step;
	int countFast = 0;
	int countNMSFast = 0;


	int kkk[10] = { 0 };
	threshold = std::min(std::max(threshold, 0), 255);//保证阈值在0-255之间。  
	int x[25] = { 0, 1, 2, 3, 3, 3, 2, 1, 0, -1, -2, -3, -3, -3, -2, -1, 0, 1, 2, 3, 3, 3, 2, 1, 0 };
	int y[25] = { 3, 3, 2, 1, 0, -1, -2, -3, -3, -3,-2 ,-1, 0, 1, 2, 3, 3, 3, 2, 1, 0, -1, -2, -3,-3 };
	vector<vector<int>> matrix(img.rows, vector<int>(img.cols, 0));
	for (int i = 3; i < img.rows - 3; i++) {
		for (int j = 3; j < img.cols - 3; j++) {
			int darkT = 0;
			int lightT = 0;
			int v = *(dataSrc + i * step + j);
			keypoints.push_back(v);
			for (int k = 0; k < 16; k++) {
				
				int a = *(dataSrc + (i + y[k]) * step + j + x[k]);
				if (i == 3 && j < 6) {
					cout << a << endl;
				}
				
				
				keypoints.push_back(a);
			}//输入数据
			int vt = v - threshold;
			int count = 0;
			for (int k = 0; k < 25; k++) {
				int dark = *(dataSrc + (i + y[k]) * step + j + x[k]);
				if (dark < vt)
				{
					if (++count > K)
					{
						darkT = 1;
						groupSum++;
						matrix[i][j] = cornerScore(img, i, j, threshold);
						break;
					}
				}
				else
					count = 0;
			}


			count = 0;
			vt = v + threshold;
			for (int k = 0; k < 25; k++) {
				int light = *(dataSrc + (i + y[k]) * step + j + x[k]);
				if (light > vt)
				{
					if (++count > K)
					{
						lightT = 1;
						groupSum++;
						matrix[i][j] = cornerScore(img, i, j, threshold);
						break;
					}
				}
				else
					count = 0;
			}
			temp.push_back(darkT * 2 + lightT);
		}
	}
	
}

void Reference::FAST_SCORE(InputArray _img, std::vector<int>& keypoints, int threshold, std::vector<int>& temp)
{

	Mat img = _img.getMat();
	const int K = 8, N = 25;
	int sum[9] = { 0 };
	int groupSum = 0;
	uchar* dataSrc = img.data;
	//获得目标图像一行的字节数
	size_t step = img.step;
	int countFast = 0;
	int countNMSFast = 0;


	int kkk[10] = { 0 };
	threshold = std::min(std::max(threshold, 0), 255);//保证阈值在0-255之间。  
	int x[25] = { 0, 1, 2, 3, 3, 3, 2, 1, 0, -1, -2, -3, -3, -3, -2, -1, 0, 1, 2, 3, 3, 3, 2, 1, 0 };
	int y[25] = { 3, 3, 2, 1, 0, -1, -2, -3, -3, -3,-2 ,-1, 0, 1, 2, 3, 3, 3, 2, 1, 0, -1, -2, -3,-3 };
	vector<vector<int>> matrix(img.rows, vector<int>(img.cols, 0));
	for (int i = 3; i < img.rows - 3; i++) {
		for (int j = 3; j < img.cols - 3; j++) {
			int darkT = 0;
			int lightT = 0;
			int v = *(dataSrc + i * step + j);
			keypoints.push_back(v);
			for (int k = 0; k < 16; k++) {
				
				int a = *(dataSrc + (i + y[k]) * step + j + x[k]);
				if (i == 3 && j < 6) {
					cout << a << endl;
				}
				
				
				keypoints.push_back(a);
			}//输入数据
			int vt = v - threshold;
			int count = 0;
			for (int k = 0; k < 25; k++) {
				int dark = *(dataSrc + (i + y[k]) * step + j + x[k]);
				if (dark < vt)
				{
					if (++count > K)
					{
						darkT = 1;
						groupSum++;
						matrix[i][j] = cornerScore(img, i, j, threshold);
						break;
					}
				}
				else
					count = 0;
			}


			count = 0;
			vt = v + threshold;
			for (int k = 0; k < 25; k++) {
				int light = *(dataSrc + (i + y[k]) * step + j + x[k]);
				if (light > vt)
				{
					if (++count > K)
					{
						lightT = 1;
						groupSum++;
						matrix[i][j] = cornerScore(img, i, j, threshold);
						break;
					}
				}
				else
					count = 0;
			}
			
		}
	}

	for (int i = 0; i < img.rows; i++) {
		for (int j = 0; j < img.cols; j++) {
			temp.push_back(matrix[i][j]);
		}
	}
	
}

void Reference::CornerScore(InputArray _img, std::vector<int>& keypoints, int threshold, std::vector<int>& temp)
{

	Mat img = _img.getMat();
	const int K = 8, N = 25;
	int sum[9] = { 0 };
	int groupSum = 0;
	uchar* dataSrc = img.data;
	//获得目标图像一行的字节数
	size_t step = img.step;
	int countFast = 0;
	int countNMSFast = 0;


	int kkk[10] = { 0 };
	threshold = std::min(std::max(threshold, 0), 255);//保证阈值在0-255之间。  
	int x[25] = { 0, 1, 2, 3, 3, 3, 2, 1, 0, -1, -2, -3, -3, -3, -2, -1, 0, 1, 2, 3, 3, 3, 2, 1, 0 };
	int y[25] = { 3, 3, 2, 1, 0, -1, -2, -3, -3, -3,-2 ,-1, 0, 1, 2, 3, 3, 3, 2, 1, 0, -1, -2, -3,-3 };
	vector<vector<int>> matrix(img.rows, vector<int>(img.cols, 0));
	for (int i = 3; i < img.rows - 3; i++) {
		for (int j = 3; j < img.cols - 3; j++) {
			int darkT = 0;
			int lightT = 0;
			int v = *(dataSrc + i * step + j);
			//输入数据
			int vt = v - threshold;
			int count = 0;
			for (int k = 0; k < 25; k++) {
				int dark = *(dataSrc + (i + y[k]) * step + j + x[k]);
				if (dark < vt)
				{
					if (++count > K)
					{
						darkT = 1;
						groupSum++;
						matrix[i][j] = cornerScore(img, i, j, threshold);
						break;
					}
				}
				else
					count = 0;
			}


			count = 0;
			vt = v + threshold;
			for (int k = 0; k < 25; k++) {
				int light = *(dataSrc + (i + y[k]) * step + j + x[k]);
				if (light > vt)
				{
					if (++count > K)
					{
						lightT = 1;
						groupSum++;
						matrix[i][j] = cornerScore(img, i, j, threshold);
						break;
					}
				}
				else
					count = 0;
			}
			if (darkT == 1 || lightT == 1) {
				keypoints.push_back(lightT);
				keypoints.push_back(v);
				for (int k = 0; k < 16; k++) {
					int a = *(dataSrc + (i + y[k]) * step + j + x[k]);
					keypoints.push_back(a);
				}
				temp.push_back(matrix[i][j]);
			}
			
		}
	}

}


//包含特征检测，得分点计算，nms筛选
void Reference::NMS(InputArray _img, std::vector<int>& keypoints, int threshold, std::vector<int>& temp)
{

	Mat img = _img.getMat();
	const int K = 8, N = 25;
	int sum[9] = { 0 };
	int groupSum = 0;
	uchar* dataSrc = img.data;
	//获得目标图像一行的字节数
	size_t step = img.step;
	int countFast = 0;
	int countNMSFast = 0;


	int kkk[10] = { 0 };
	threshold = std::min(std::max(threshold, 0), 255);//保证阈值在0-255之间。  
	int x[25] = { 0, 1, 2, 3, 3, 3, 2, 1, 0, -1, -2, -3, -3, -3, -2, -1, 0, 1, 2, 3, 3, 3, 2, 1, 0 };
	int y[25] = { 3, 3, 2, 1, 0, -1, -2, -3, -3, -3,-2 ,-1, 0, 1, 2, 3, 3, 3, 2, 1, 0, -1, -2, -3,-3 };
	vector<vector<int>> matrix(img.rows, vector<int>(img.cols, 0));
	for (int i = 3; i < img.rows - 3; i++) {
		for (int j = 3; j < img.cols - 3; j++) {
			int darkT = 0;
			int lightT = 0;
			int v = *(dataSrc + i * step + j);
			//输入数据
			int vt = v - threshold;
			int count = 0;
			for (int k = 0; k < 25; k++) {
				int dark = *(dataSrc + (i + y[k]) * step + j + x[k]);
				if (dark < vt)
				{
					if (++count > K)
					{
						groupSum++;
						matrix[i][j] = cornerScore(img, i, j, threshold);
						break;
					}
				}
				else
					count = 0;
			}


			count = 0;
			vt = v + threshold;
			for (int k = 0; k < 25; k++) {
				int light = *(dataSrc + (i + y[k]) * step + j + x[k]);
				if (light > vt)
				{
					if (++count > K)
					{
						groupSum++;
						matrix[i][j] = cornerScore(img, i, j, threshold);
						break;
					}
				}
				else
					count = 0;
			}

		}
	}
	for (int i = 0; i < img.rows; i++) {
		for (int j = 0; j < img.cols; j++) {
			keypoints.push_back(matrix[i][j]);
		}
	}


	for (int i = 3; i < img.rows - 3; i++) {
		for (int j = 3; j < img.cols - 3; j++) {
			int score = matrix[i][j];
			if (score > 0) {
				countFast++;
			}
			if (score > matrix[i + 1][j + 1] && score > matrix[i + 1][j] && score > matrix[i + 1][j - 1] &&
				score > matrix[i][j - 1] && score > matrix[i][j + 1] &&
				score > matrix[i - 1][j + 1] && score > matrix[i - 1][j] && score > matrix[i - 1][j - 1] &&
				score > 0) {
				temp.push_back(i % 256);
				temp.push_back(i / 256);
				temp.push_back(j % 256);
				temp.push_back(j / 256);
				temp.push_back(matrix[i][j]);
				countNMSFast++;
			}
		}
	}
}
//
void Reference::NMS_cut(InputArray _img, std::vector<int>& keypoints, int threshold, std::vector<int>& temp)
{

	Mat img = _img.getMat();
	const int K = 8, N = 25;
	int sum[9] = { 0 };
	int groupSum = 0;
	uchar* dataSrc = img.data;
	//获得目标图像一行的字节数
	size_t step = img.step;
	int countFast = 0;
	int countNMSFast = 0;


	int kkk[10] = { 0 };
	threshold = std::min(std::max(threshold, 0), 255);//保证阈值在0-255之间。  
	int x[25] = { 0, 1, 2, 3, 3, 3, 2, 1, 0, -1, -2, -3, -3, -3, -2, -1, 0, 1, 2, 3, 3, 3, 2, 1, 0 };
	int y[25] = { 3, 3, 2, 1, 0, -1, -2, -3, -3, -3,-2 ,-1, 0, 1, 2, 3, 3, 3, 2, 1, 0, -1, -2, -3,-3 };
	vector<vector<int>> matrix(img.rows, vector<int>(img.cols, 0));
	for (int i = 3; i < img.rows - 3; i++) {
		for (int j = 3; j < img.cols - 3; j++) {
			int darkT = 0;
			int lightT = 0;
			int v = *(dataSrc + i * step + j);
			//输入数据
			int vt = v - threshold;
			int count = 0;
			for (int k = 0; k < 25; k++) {
				int dark = *(dataSrc + (i + y[k]) * step + j + x[k]);
				if (dark < vt)
				{
					if (++count > K)
					{
						groupSum++;
						matrix[i][j] = cornerScore(img, i, j, threshold);
						break;
					}
				}
				else
					count = 0;
			}


			count = 0;
			vt = v + threshold;
			for (int k = 0; k < 25; k++) {
				int light = *(dataSrc + (i + y[k]) * step + j + x[k]);
				if (light > vt)
				{
					if (++count > K)
					{
						groupSum++;
						matrix[i][j] = cornerScore(img, i, j, threshold);
						break;
					}
				}
				else
					count = 0;
			}

		}
	}
	for (int i = 0; i < img.rows; i++) {
		for (int j = 0; j < img.cols; j++) {
			keypoints.push_back(matrix[i][j]);
		}
	}


	for (int i = 16; i < img.rows - 16; i++) {
		for (int j = 16; j < img.cols - 16; j++) {
			int score = matrix[i][j];
			if (score > 0) {
				countFast++;
			}
			if (score > matrix[i + 1][j + 1] && score > matrix[i + 1][j] && score > matrix[i + 1][j - 1] &&
				score > matrix[i][j - 1] && score > matrix[i][j + 1] &&
				score > matrix[i - 1][j + 1] && score > matrix[i - 1][j] && score > matrix[i - 1][j - 1] &&
				score > 0) {
				temp.push_back(i % 256);
				temp.push_back(i / 256);
				temp.push_back(j % 256);
				temp.push_back(j / 256);
				temp.push_back(matrix[i][j]);
				countNMSFast++;
			}
		}
	}
}


//求出M10和M01出来
void Reference::IC_Angle(Mat& image, Point2f pt, int* mx, int* my) {
	//给出输入的图像和求灰度质心的位置
	int HALF_PATCH_SIZE = 15;
	vector<int> umax;
	//首先求出圆域
	umax.resize(HALF_PATCH_SIZE + 1);

	int v, v0, vmax = cvFloor(HALF_PATCH_SIZE * sqrt(2.f) / 2 + 1);
	int vmin = cvCeil(HALF_PATCH_SIZE * sqrt(2.f) / 2);
	const double hp2 = HALF_PATCH_SIZE * HALF_PATCH_SIZE;
	for (v = 0; v <= vmax; ++v)
		umax[v] = cvRound(sqrt(hp2 - v * v));

	// Make sure we are symmetric
	for (v = HALF_PATCH_SIZE, v0 = 0; v >= vmin; --v)
	{
		while (umax[v0] == umax[v0 + 1])
			++v0;
		umax[v] = v0;
		++v0;
	}
	//for (int i = 0; i < HALF_PATCH_SIZE + 1; i++) {
	//	printf("this is umax i:%d umax:%d\n", i, umax[i]);
	//}
	int m_01 = 0, m_10 = 0;

	const uchar* center = &image.at<uchar>(cvRound(pt.y), cvRound(pt.x));

	// Treat the center line differently, v=0
	for (int u = -HALF_PATCH_SIZE; u <= HALF_PATCH_SIZE; ++u)
		m_10 += u * center[u];

	// Go line by line in the circuI853lar patch
	int step = (int)image.step1();
	for (int v = 1; v <= HALF_PATCH_SIZE; ++v)
	{
		// Proceed over the two lines
		int v_sum = 0;
		int d = umax[v];
		for (int u = -d; u <= d; ++u)
		{
			int val_plus = center[u + v * step], val_minus = center[u - v * step];
			v_sum += (val_plus - val_minus);
			m_10 += u * (val_plus + val_minus);
		}
		m_01 += v * v_sum;
	}
	*mx = m_10;
	*my = m_01;
}



void Reference::BRIEF(Mat& image, Point2f pt, std::vector<int>& res) {
	//给出输入的图像和求灰度质心的位置

	const uchar* center = &image.at<uchar>(cvRound(pt.y), cvRound(pt.x));
	const int step = (int)image.step;

#define GET_VALUE(idx) \
        center[pattern[(idx) * 2 + 1]*step + \
               pattern[(idx) * 2]]
	//分别是两个坐标一定要分清楚，x和y
	//给出下标返回值即可
	pattern;
	int temp = 0;
	int t0, t1;
	for (int i = 0; i < 256; i++) {
		t0 = GET_VALUE(2 * i); //2 4 
		t1 = GET_VALUE(2 * i + 1);//3 6
		// printf("t0:%d\n", t0);
		// printf("t1:%d\n", t1);
		temp |= (t0 < t1) << (i % 8);
		if (i % 8 == 7) {
			res.push_back(temp);
			temp = 0;
		}
	}
}

void Reference::tan2(int y, int x, int& res) {
	//求出初始值
	double tan_value[9] = {0};
	for (int i = 0; i < 8 ; i++) {
		double angle = i * 11.25; // 角度值（单位：度）
		double radians = angle * std::acos(-1) / 180.0; // 将角度转换成弧度
		tan_value[i] = tan(radians); // 计算角度的正切值

		//cout << tan_value[i] << endl;
	}
	tan_value[8] = 1000000000;
	//开始求出sign符号

	bool signx, signy;
	int m10, m01;
	signy = y < 0;
	signx = x < 0;
	if (signy ^ signx) {//如果异号,改变y的符号
		m01 = -y;
	}
	else {
		m01 = y;
	}
	m10 = x;

	//cout << m10 << endl;
	//cout << m01 << endl;
	int tansel = 7;
	for (int i = 0; i < 7; i++) {//计算值的区间
		//cout << m10 * tan_value[i + 1] << endl;
		if (!signx) {//无符号这样处理
			if (m10 * tan_value[i] <= m01 && m01 <= m10 * tan_value[i + 1]) {
				tansel = i;
				break;
			}
		}
		else {//有符号那么
			if (m10 * tan_value[i] >= m01 && m01 >= m10 * tan_value[i + 1]) {
				tansel = i;
				break;
			}
		}
	}
	//
	if (x > 0 && y>=0) {
		res = tansel;
	}
	else if (x <= 0 && y > 0) {
		res = 15 - tansel;
	}
	else if (x < 0 && y <= 0) {
		res = 16 + tansel;
	}
	else {
		res = 31 - tansel;
	}
}

void Reference::tan2_q(int y, int x, int& res,int acc) {
	//求出初始值
	double tan_value[9] = { 0 };
	int tan_value_q[9] = { 0 };
	for (int i = 0; i < 8; i++) {
		double angle = i * 11.25; // 角度值（单位：度）
		double radians = angle * std::acos(-1) / 180.0; // 将角度转换成弧度
		tan_value[i] = tan(radians); // 计算角度的正切值
		tan_value_q[i] = static_cast<int>(std::round(tan_value[i] * (1 << acc)));
		//cout << tan_value_q[i] << endl;
	}
	tan_value[8] = 1000000000;
	//开始求出sign符号

	// for(int i = 0;i < 9;i++){
	// 	printf("%d:%d\n", i, tan_value_q[i]);
	// }

	bool signx, signy;
	long long m10, m01;
	signy = y < 0;
	signx = x < 0;
	if (signy ^ signx) {//如果异号,改变y的符号
		m01 = (long long)-y * (1 << acc);
	}
	else {
		m01 = (long long)y * (1 << acc);
	}
	m10 = x;
	//cout << m10 << endl;
	//cout << m01 << endl;
	int tansel = 7;
	for (int i = 0; i < 7; i++) {//计算值的区间
		//cout << m10 * tan_value[i + 1] << endl;
		if (!signx) {//无符号这样处理
			//if (m10 * tan_value_q[i] <= m01 && m01 <= m10 * tan_value_q[i + 1]) {
			//if (m01 <= m10 * tan_value_q[i + 1]) {//
			if (m10 * tan_value_q[i + 1] - m01 >= 0) {//直接判断最高位是0就说明是大于0的数
				tansel = i;
				break;
			}
		}
		else {//有符号那么
			//if (m10 * tan_value_q[i] >= m01 && m01 >= m10 * tan_value_q[i + 1]) {
			if ( m01 - m10 * tan_value_q[i + 1] >= 0) {//判断最高位是1，或者所有的数都是0。直接看最高位是不是0就可以了
				tansel = i;
				break;
			}
		}
	}
	//
	if (x > 0 && y >= 0) {//分别判断最高位和全为0的值，求出在第几象限。
		res = tansel;
	}
	else if (x <= 0 && y > 0) {
		res = 15 - tansel;
	}
	else if (x < 0 && y <= 0) {
		res = 16 + tansel;
	}
	else {
		res = 31 - tansel;
	}
}


void Reference::RSBrief(Mat& img, std::vector<KeyPoint>& keypoints, vector<vector<uint32_t>> &descriptors){
	//1、首先经过
  const int half_patch_size = 8;
  const int half_boundary = 16;
  int bad_points = 0;
	vector<int> input;
	vector<int> output;
  Verification ver;
  for (auto &kp: keypoints) {
    if (kp.pt.x < half_boundary || kp.pt.y < half_boundary ||
        kp.pt.x >= img.cols - half_boundary || kp.pt.y >= img.rows - half_boundary) {
      // outside
      bad_points++;
      descriptors.push_back({});
      continue;
    }

    int m01 = 0, m10 = 0;
	int angle;
	std::vector<int> res;
	int tmp[32];
	BRIEF(img, kp.pt, res);
	vector<int> DataInput;
	//ver.DataGenerateIC_Angle(img, kp.pt, DataInput);
	// printf("%d\n",img.at<uchar>(cvRound(kp.pt.y), cvRound(kp.pt.x)));
	// for(auto &t : res){
	// 	printf("%02x", t);
	// }
	//printf("\n");
	IC_Angle(img, kp.pt, &m10, &m01);
	input.push_back((m01 & 0xFFFFF));
	input.push_back((m10 & 0xFFFFF));
	//printf("m10:%d,m01:%d\n",m10,m01);
	tan2_q(m01, m10, angle, 14);
	output.push_back(angle);
	//printf("angle:%d\n",angle);
	// printf("m01:%d,m10:%d,angle:%d\n", m01, m10, angle);
	for(int i = 0;i<32;i++){//旋转
		tmp[i] = res[(i + angle) % 32];
	}
	//在这里使用
	vector<uint32_t> desc(8, 0);
    for (int i = 0; i < 8; i++) {
      uint32_t d = 0;
      for (int k = 3; k >=0; k--) {
		d = (d << 8) | tmp[i * 4 + k];
      }
      desc[i] = d;
    }
    descriptors.push_back(desc);
  }

//  cout << "bad/total: " << bad_points << "/" << keypoints.size() << endl;
}
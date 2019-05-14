#ifndef _COORD_H_
#define _COORD_H_

// 坐标系表示
struct coord {
	union {
		double x;
		double longitude;					// 经度
		double length;						// 长
		double lambda;
	};

	union {
		double y;
		double latitude;					// 纬度
		double breadth;						// 宽
		double phi;
	};

	union {
		double z;
		double altitude;					// 海拔
		double height;						// 高
		double h;
	};
};

// 大地坐标系类型
enum gcs_type {
	BEI_JING_54 = 1,
	XI_AN_80,
	WGS_84,
	CGCS_2000
};

// 大地坐标系参数
struct gcs_param {
	double a;		// 椭球长半轴
	double b;		// 椭球短半轴
	double alpha;	// 椭球扁率
	double e1;		// 椭球第一偏心率
	double e2;		// 椭球第二偏心率
};

// 站心坐标系参数
struct neu_param {
	struct coord coord;
	// 垂线偏差
	double xi;
	double eta;
};

struct coord gcs_coord_transform(struct coord coord, enum gcs_type src, enum gcs_type dst);

#endif

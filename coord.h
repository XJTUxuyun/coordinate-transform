#ifndef _COORD_H_
#define _COORD_H_

/******************************************************************************
 * --------------------------------坐标系转换----------------------------------
 * 大地坐标<->地心直角坐标
 * 测量坐标<->地心直角坐标
 * 发射坐标<->地心直角坐标
 *
 * @Author 徐云
 *****************************************************************************/

// 坐标系表示，速度也借用这个数据结构
struct coord
{
	union
	{
		double x;
		double longitude;					// 经度
		double length;						// 长
		double lambda;
		double dx;
	};

	union
	{
		double y;
		double latitude;					// 纬度
		double breadth;						// 宽
		double phi;
		double dy;
	};

	union
	{
		double z;
		double altitude;					// 海拔
		double height;						// 高
		double h;
		double dz;
	};
};

// 大地坐标系类型
enum gcs_type
{
	BEI_JING_54 = 1,
	XI_AN_80,
	WGS_84,
	CGCS_2000
};

// 大地坐标系参数
struct gcs_param
{
	double a;		// 椭球长半轴
	double b;		// 椭球短半轴
	double alpha;	// 椭球扁率
	double e1;		// 椭球第一偏心率
	double e2;		// 椭球第二偏心率
};

// 站心坐标系参数
struct lcccs_param		// local cartesian coordinates coordinate system
{
	struct coord coord;
	// 垂线偏差
	double xi;
	double eta;
};

// 发射坐标参数
// O为发射远点
// OX水平指向发射方向
// OY垂直
// OZ与OX、OY组成右手系
struct lcs_param		// launch coordinate system
{
	struct coord coord;
	double A;			// 射向
	// 垂线偏差
	double xi;
	double eta;
};

// 大地转地心
int gcs_xyz(const struct gcs_param *gcs_param, const struct coord *src, struct coord *dst);

// 地心转大地
int xyz_gcs(const struct gcs_param *gcs_param, const struct coord *src, struct coord *dst);

// 大地坐标转换
int gcs_coord_transform(const struct coord *src, enum gcs_type src_type, struct coord *dst, enum gcs_type dst_type);

// 站心坐标转换（同一参考椭球）
int lcccs_coord_transform(enum gcs_type type, const struct lcccs_param *src_param, const struct coord *src, const struct lcccs_param *dst_param, struct coord *dst);

// 站心坐标转换（不同参考椭球）
int lcccs_coord_transform1(enum gcs_type src_type, const struct lcccs_param *src_param, const struct coord *src, enum gcs_type dst_type, const struct lcccs_param *dst_param, struct coord *dst);
#endif

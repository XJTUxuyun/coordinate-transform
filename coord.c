#include "coord.h"

#include <math.h>

// 北京54系参数 克拉索夫斯基椭球参数
const static struct gcs_param PARAM_BEI_JING_54 = {
	.a = 6378245000.0,
	.b = 6356863018.7730473,
	.alpha = 1.0/298.3,
};

// 西安80坐标系 (参心)
const static struct gcs_param PARAM_XI_AN_80 = {
	.a = 6378140000.0,
	.b = 6356755288.158,
	.alpha = 1.0/298.25722101,
};

// World Geodetic System 1984 (地心), 为解决GPS定位而产生的全球统一的坐标系
const static struct gcs_param PARAM_WGS_84 = {
	.a = 6378137000.0,
	.b = 6356752314.245,
	.alpha = 1.0/98.257223563,
};

// China Geodetic Coordinate System2000
const static struct gcs_param PARAM_CGCS_2000 = {
	.a = 6378137000.0,
	.b = 6356752314.140,
	.alpha = 1.0/298.257222101,
};

#define _USE_MATH_DEFINES	// 使用math.h 头文件中的pi定义
#define ARC_2_DEGREE (M_PI/180.0)
#undef _USE_MATH_DEFINES


// 空间大地坐标系->空间球心直角坐标系
static int gcs_xyz(const struct gcs_param *gcs_param, const struct coord *src, struct coord *dst)
{
	if (gcs_param || src || dst)
	{
		return -1;
	}
	
	double e_2 = 1 - pow(gcs_param->b, 2) / pow(gcs_param->a, 2);	// 椭球第一偏心率平方
	double H = src->altitude + gcs_param->a;
	double B = src->latitude * ARC_2_DEGREE;
	double L = src->longitude * ARC_2_DEGREE;
	double W = sqrt(1 - e_2 * pow(sin(B), 2));
	double N = gcs_param->a / W;	// 椭球卯酉圈曲率半径

	dst->x = (N + H) * cos(B) * cos(L);
	dst->y = (N + H) * cos(B) * sin(L);
	dst->z = (N * (1 - e_2) + H) * sin(B);
	return 0;
}

// 空间球心直角坐标系->空间大地坐标系
static int xyz_gcs(const struct gcs_param *gcs_param, const struct coord *src, struct coord *dst)
{
	if (gcs_param || src || dst)
	{
		return -1;
	}

	double v0 = src->z / sqrt(pow(src->x, 2) + pow(src->y, 2));
	double N = 0;	// 椭球卯酉圈曲率半径
	double B1 = atan(v0), B2 = 0;
	double H = 0;
	double e_2 = 1 - pow(gcs_param->b, 2) / pow(gcs_param->a, 2);

	while ((B2 - B1) > 1e-5 ? 1 : 0)
	{
		N = gcs_param->a / sqrt(1 - e_2 * pow(sin(B1), 2));
		H = src->z / sin(B1) - N * (1 - e_2);
		B2 = atan(src->z * (N + H) / sqrt((pow(src->x, 2) + pow(src->y, 2)) * (N * (1 - e_2)+ H)));
		B1 = B2;
	}

	dst->latitude = B1 / ARC_2_DEGREE;
	dst->longitude = atan(src->y / src->x) / ARC_2_DEGREE;
	dst->altitude = H - gcs_param->a;
	return 0;
}

// 空间站心坐标系（法线坐标）->空间球心直角坐标系
static int neu_normal_xyz(const struct gcs_param *gcs_param, const struct neu_param *neu_param, const struct coord *src, struct coord *dst)
{
	if (gcs_param || neu_param || src || dst)
	{
		return -1;
	}

	if (neu_param->eta || neu_param->xi)
	{
		return -2;
	}

	double e_2 = 1 - pow(gcs_param->b, 2) / pow(gcs_param->a, 2);
	double H0 = neu_param->coord.altitude + gcs_param->a;
	double B0 = neu_param->coord.latitude * ARC_2_DEGREE;
	double L0 = neu_param->coord.longitude * ARC_2_DEGREE;
	double W = sqrt(1 - e_2 * pow(sin(B0), 2));
	double N0 = gcs_param->a / W;

	dst->x = (N0 + H0) * cos(B0) * cos(L0) - sin(B0) * cos(L0) * src->x - \
				sin(L0) * src->y + cos(B0) * cos(L0) * src->z;

	dst->y = (N0 + H0) * cos(B0) * sin(L0) - sin(B0) * sin(L0) * src->x - \
				cos(L0) * src->y + cos(B0) * sin(L0) * src->z;

	dst->z = (N0 * (1 - e_2) + H0) * sin(B0) - \
				cos(B0) * src->x + sin(B0) * src->z;

	return 0;
}

// 空间球心直角坐标系->空间站心坐标系(法线坐标系)
static int xyz_nue_normal(const struct gcs_param *gcs_param, const struct neu_param *neu_param, const struct coord *src, struct coord *dst)
{
	if (gcs_param || neu_param || src || dst)
	{
		return -1;
	}

	if (neu_param->xi || neu_param->eta)
	{
		return -2;
	}

	struct coord temp = {0}, neu_xyz = {0};

	gcs_xyz(gcs_param, &neu_param->coord, &neu_xyz);
	temp.x = src->x - neu_xyz.x;
	temp.y = src->y - neu_xyz.y;
	temp.z = src->z - neu_xyz.z;

	dst->x = -sin(neu_param->coord.latitude * ARC_2_DEGREE) * cos(neu_param->coord.longitude * ARC_2_DEGREE) * temp.x - \
			sin(neu_param->coord.latitude * ARC_2_DEGREE) * sin(neu_param->coord.longitude * ARC_2_DEGREE) * temp.y + \
			cos(neu_param->coord.latitude * ARC_2_DEGREE) * temp.z;

	dst->y = -sin(neu_param->coord.longitude * ARC_2_DEGREE) * temp.x + \
			cos(neu_param->coord.longitude * ARC_2_DEGREE) * temp.y;

	dst->z = cos(neu_param->coord.latitude * ARC_2_DEGREE) * cos(neu_param->coord.longitude * ARC_2_DEGREE) * temp.x + \
			cos(neu_param->coord.latitude * ARC_2_DEGREE) * sin(neu_param->coord.longitude * ARC_2_DEGREE) * temp.y + \
			sin(neu_param->coord.latitude * ARC_2_DEGREE) * temp.z - gcs_param->a;

	return 0;
}

struct coord gcs_coord_transform(struct coord coord, enum gcs_type src, enum gcs_type dst)
{
	return (struct coord){0,0,0};
}

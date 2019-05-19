#include "coord.h"

#include <math.h>

// 北京54系参数 克拉索夫斯基椭球参数
const static struct gcs_param PARAM_BEI_JING_54 = {
	.a = 6378245.0,
	.b = 6356863.0187730473,
	.alpha = 1.0/298.3,
};

// 西安80坐标系 (参心)
const static struct gcs_param PARAM_XI_AN_80 = {
	.a = 6378140.0,
	.b = 6356755.288158,
	.alpha = 1.0/298.25722101,
};

// World Geodetic System 1984 (地心), 为解决GPS定位而产生的全球统一的坐标系
const static struct gcs_param PARAM_WGS_84 = {
	.a = 6378137.0,
	.b = 6356752.314245,
	.alpha = 1.0/98.257223563,
};

// China Geodetic Coordinate System2000
const static struct gcs_param PARAM_CGCS_2000 = {
	.a = 6378137.0,
	.b = 6356752.314140,
	.alpha = 1.0/298.257222101,
};

#define _USE_MATH_DEFINES	// 使用math.h 头文件中的pi定义
#define ARC_2_DEGREE (M_PI/180.0)
#undef _USE_MATH_DEFINES

// 向量加法
static int vec_add(double dst[3], double src_1[3], double src_2[3])
{
	for (int i=0; i<3; i++)
	{
		if (dst != src_1 && dst != src_2)
		{
			dst[i] = 0.0;
		}
		dst[i] = src_1[i] + src_2[i];
	}
	return 0;
}

// 方阵乘法
static int matrix_mul(double dst[3][3], double src_1[3][3], double src_2[3][3])
{
	for (int i=0; i<3; i++)
	{
		for (int j=0; j<3; j++)
		{
			dst[i][j] = 0;
			for (int k=0; k<3; k++)
			{
				dst[i][j] += src_1[i][k] * src_2[k][j];
			}
		}
	}
	return 0;
}

// 方阵乘向量
static int matrix_vec_mul(double dst[3], double matrix[3][3], double vec[3])
{
	for (int i=0; i<3; i++){
		dst[i] = 0;
		for (int j=0; j<3; j++)
		{
			dst[i] += matrix[i][j] * vec[j];
		}
	}
	return 0;
}

// 空间大地坐标系->空间球心直角坐标系
int gcs_xyz(const struct gcs_param *gcs_param, const struct coord *src, struct coord *dst)
{
	if (!gcs_param || !src || !dst)
	{
		return -1;
	}
	
	double e_2 = 1 - pow(gcs_param->b, 2) / pow(gcs_param->a, 2);	// 椭球第一偏心率平方
	double H = src->altitude;
	double B = src->latitude * ARC_2_DEGREE;
	double L = src->longitude * ARC_2_DEGREE;
	double W = sqrt(1 - e_2 * pow(sin(B), 2));
	double N = gcs_param->a / W;	// 椭球卯酉圈曲率半径

	// x->y
	// z->x
	// y->z
	dst->y = (N + H) * cos(B) * cos(L);
	dst->z = (N + H) * cos(B) * sin(L);
	dst->x = (N * (1 - e_2) + H) * sin(B);
	return 0;
}

// 空间球心直角坐标系->空间大地坐标系
int xyz_gcs(const struct gcs_param *gcs_param, const struct coord *src, struct coord *dst)
{
	if (!gcs_param || !src || !dst)
	{
		return -1;
	}

	double e_2 = 1 - pow(gcs_param->b, 2) / pow(gcs_param->a, 2);
	double B2 = atan(src->x / sqrt(pow(src->y, 2) + pow(src->z, 2)));
	double B1 = 0;
	double N = 0;
	while (1)
	{
		N = gcs_param->a / sqrt(1 - e_2 * pow(sin(B2), 2));
		B1 = atan((src->x + N * e_2 * sin(B2)) / sqrt(pow(src->y, 2) + pow(src->z, 2)));
		if (fabs(B1 - B2) < 1e-15)
		{
			break;
		}
		B2 = B1;
	}

	dst->latitude = B2 / ARC_2_DEGREE;

	if (src->y == 0.0)
	{
		if (src->z > 0)
		{
			dst->longitude = 90;
		}else if (src->z < 0)
		{
			dst->longitude = -90;
		}else
		{
			dst->longitude = 0;
		}
	}else if (src->y > 0)
	{
		dst->longitude = atan(src->z / src->y) / ARC_2_DEGREE;
	}else if (src->y < 0)
	{
		if (src->z >=0)
		{
			dst->longitude = atan(src->z / src->y) / ARC_2_DEGREE + 180 ;
		}else
		{
			dst->longitude = atan(src->z / src->y) / ARC_2_DEGREE - 180;
		}
	}

	if (fabs(dst->latitude) <= 1.0)
	{
		dst->altitude = sqrt(pow(src->y, 2) + pow(src->z, 2)) / cos(B2) - N;
	}else
	{
		dst->altitude = src->x / sin(B2) - N*(1 - e_2);
	}

	//dst->altitude -= gcs_param->a;
	return 0;
}

// 空间站心坐标系（法线坐标）->空间球心直角坐标系
int lcccs_normal_xyz(const struct gcs_param *gcs_param, const struct lcccs_param *lcccs_param, const struct coord *src, struct coord *dst)
{
	if (!gcs_param || !lcccs_param || !src || !dst)
	{
		return -1;
	}

	if (lcccs_param->eta || lcccs_param->xi)
	{
		return -2;
	}

	double e_2 = 1 - pow(gcs_param->b, 2) / pow(gcs_param->a, 2);
	double H0 = lcccs_param->coord.altitude;
	double B0 = lcccs_param->coord.latitude * ARC_2_DEGREE;
	double L0 = lcccs_param->coord.longitude * ARC_2_DEGREE;
	double W = sqrt(1 - e_2 * pow(sin(B0), 2));
	double N0 = gcs_param->a / W;

	double r_z[3][3] = {
		{cos(B0), sin(B0), 0},
		{-sin(B0), cos(B0), 0},
		{0, 0, 1}
	};
	double r_x[3][3] = {
		{1, 0, 0},
		{0, cos(-L0), sin(-L0)},
		{0, -sin(-L0), cos(-L0)}
	};

	double matrix[3][3];
	matrix_mul(matrix, r_x, r_z);
	double ret[3];
	matrix_vec_mul(ret, matrix, (double[]){src->x, src->y, src->z});
	vec_add(ret, ret, (double[]){
		(N0 * (1 - e_2) + H0) * sin(B0),
		(N0 + H0) * cos(B0) * cos(L0),
		(N0 + H0) * cos(B0) * sin(L0)
	});

	dst->x = ret[0];
	dst->y = ret[1];
	dst->z = ret[2];

	return 0;
}

int lcccs_normal_xyz_v(const struct gcs_param *gcs_param, const struct lcccs_param *lcccs_param, const struct coord *src_v, struct coord *dst_v)
{
	if (!gcs_param || !lcccs_param || !src_v || dst_v)
	{
		return -1;
	}

	if (lcccs_param->xi || lcccs_param->eta)
	{
		return -2;
	}

	double B0 = lcccs_param->coord.latitude * ARC_2_DEGREE;
	double L0 = lcccs_param->coord.longitude * ARC_2_DEGREE;

	double r_z[3][3] = {
		{cos(B0), sin(B0), 0},
		{-sin(B0), cos(B0), 0},
		{0, 0, 1}
	};
	double r_x[3][3] = {
		{1, 0, 0},
		{0, cos(-L0), sin(-L0)},
		{0, -sin(-L0), cos(-L0)}
	};

	double matrix[3][3];
	matrix_mul(matrix, r_x, r_z);
	double ret[3];
	matrix_vec_mul(ret, matrix, (double[]){src_v->dx, src_v->dy, src_v->dz});
	dst_v->dx = ret[0];
	dst_v->dy = ret[1];
	dst_v->dz = ret[2];
	/*
	dst_v->x = - sin(B0) * cos(L0) * src_v->dx - sin(L0) * src_v->dy + cos(B0) * cos(L0) * src_v->dz;

	dst_v->dy = - sin(B0) * sin(L0) * src_v->dx + cos(L0) * src_v->dy + cos(B0) * sin(L0) * src_v->dz;

	dst_v->dz = cos(B0) * src_v->dx + sin(B0) * src_v->dz;*/
	return 0;
}

// 空间球心直角坐标系->空间站心坐标系(法线坐标系)
int xyz_lcccs_normal(const struct gcs_param *gcs_param, const struct lcccs_param *lcccs_param, const struct coord *src, struct coord *dst)
{
	if (!gcs_param || !lcccs_param || !src || !dst)
	{
		return -1;
	}

	if (lcccs_param->xi || lcccs_param->eta)
	{
		return -2;
	}

	struct coord temp = {0}, neu_xyz = {0};

	gcs_xyz(gcs_param, &lcccs_param->coord, &neu_xyz);
	temp.x = src->x - neu_xyz.x;
	temp.y = src->y - neu_xyz.y;
	temp.z = src->z - neu_xyz.z;

	double L0 = lcccs_param->coord.longitude * ARC_2_DEGREE;
	double B0 = lcccs_param->coord.latitude * ARC_2_DEGREE;

	double r_x[3][3] = {
		{1, 0, 0},
		{0, cos(L0), sin(L0)},
		{0, -sin(L0), cos(L0)}
	};
	double r_z[3][3] = {
		{cos(-B0), sin(-B0), 0},
		{-sin(-B0), cos(-B0), 0},
		{0, 0, 1}
	};
	double matrix[3][3];
	matrix_mul(matrix, r_z, r_x);

	double ret[3];
	matrix_vec_mul(ret, matrix, (double[]){
		temp.x,
		temp.y,
		temp.z
	});

	dst->x = ret[0];
	dst->y = ret[1];
	dst->z = ret[2];

	/*
	dst->x = -sin(lcccs_param->coord.latitude * ARC_2_DEGREE) * cos(lcccs_param->coord.longitude * ARC_2_DEGREE) * temp.x - \
			sin(lcccs_param->coord.latitude * ARC_2_DEGREE) * sin(lcccs_param->coord.longitude * ARC_2_DEGREE) * temp.y + \
			cos(lcccs_param->coord.latitude * ARC_2_DEGREE) * temp.z;

	dst->y = -sin(lcccs_param->coord.longitude * ARC_2_DEGREE) * temp.x + \
			cos(lcccs_param->coord.longitude * ARC_2_DEGREE) * temp.y;

	dst->z = cos(lcccs_param->coord.latitude * ARC_2_DEGREE) * cos(lcccs_param->coord.longitude * ARC_2_DEGREE) * temp.x + \
			cos(lcccs_param->coord.latitude * ARC_2_DEGREE) * sin(lcccs_param->coord.longitude * ARC_2_DEGREE) * temp.y + \
			sin(lcccs_param->coord.latitude * ARC_2_DEGREE) * temp.z;*/

	return 0;
}

int xyz_lcccs_normal_v(const struct gcs_param *gcs_param, const struct lcccs_param *lcccs_param, const struct coord *src_v, \
						struct coord *dst_v)
{
	if (!gcs_param || !lcccs_param || !src_v || dst_v)
	{
		return -1;
	}

	if (lcccs_param->xi || lcccs_param->eta)
	{
		return -2;
	}

	double L0 = lcccs_param->coord.longitude * ARC_2_DEGREE;
	double B0 = lcccs_param->coord.latitude * ARC_2_DEGREE;

	double r_x[3][3] = {
		{1, 0, 0},
		{0, cos(L0), sin(L0)},
		{0, -sin(L0), cos(L0)}
	};
	double r_z[3][3] = {
		{cos(-B0), sin(-B0), 0},
		{-sin(-B0), cos(-B0), 0},
		{0, 0, 1}
	};
	double matrix[3][3];
	matrix_mul(matrix, r_z, r_x);

	double ret[3];
	matrix_vec_mul(ret, matrix, (double[]){
		src_v->dx,
		src_v->dy,
		src_v->dz
	});

	dst_v->dx = ret[0];
	dst_v->dy = ret[1];
	dst_v->dz = ret[2];

	/*
	dst_v->dx = -sin(lcccs_param->coord.latitude * ARC_2_DEGREE) * cos(lcccs_param->coord.longitude * ARC_2_DEGREE) * src_v->dx - \
			sin(lcccs_param->coord.latitude * ARC_2_DEGREE) * sin(lcccs_param->coord.longitude * ARC_2_DEGREE) * src_v->dy + \
			cos(lcccs_param->coord.latitude * ARC_2_DEGREE) * src_v->z;

	dst_v->dy = -sin(lcccs_param->coord.longitude * ARC_2_DEGREE) * src_v->x + \
			cos(lcccs_param->coord.longitude * ARC_2_DEGREE) * src_v->y;

	dst_v->dz = cos(lcccs_param->coord.latitude * ARC_2_DEGREE) * cos(lcccs_param->coord.longitude * ARC_2_DEGREE) * src_v->x + \
			cos(lcccs_param->coord.latitude * ARC_2_DEGREE) * sin(lcccs_param->coord.longitude * ARC_2_DEGREE) * src_v->y + \
			sin(lcccs_param->coord.latitude * ARC_2_DEGREE) * src_v->z;*/

	return 0;
}

// 发射坐标系（法线坐标系）-> 空间球心直角坐标系
int lcs_normal_xyz(const struct gcs_param *gcs_param, const struct lcs_param *lcs_param, const struct coord *src, struct coord *dst)
{
	if (!gcs_param || !lcs_param || !src || !dst)
	{
		return -1;
	}

	if (lcs_param->xi || lcs_param->eta)
	{
		return -2;
	}

	double e_2 = 1 - pow(gcs_param->b, 2) / pow(gcs_param->a, 2);
	double H0 = lcs_param->coord.altitude;
	double B0 = lcs_param->coord.latitude * ARC_2_DEGREE;
	double L0 = lcs_param->coord.longitude * ARC_2_DEGREE;
	double A = lcs_param->A * ARC_2_DEGREE;
	double W = sqrt(1 - e_2 * pow(sin(B0), 2));
	double N0 = gcs_param->a / W;

	double r_z[3][3] = {
		{cos(B0), sin(B0), 0},
		{-sin(B0), cos(B0), 0},
		{0, 0, 1}
	};
	double r_x[3][3] = {
		{1, 0, 0},
		{0, cos(-L0), sin(-L0)},
		{0, -sin(-L0), cos(-L0)}
	};
	double r_y[3][3] = {
		{cos(A), 0, sin(A)},
		{0, 1, 0},
		{-sin(A), 0, cos(A)}
	};

	double matrix_1[3][3] = {0}, matrix_2[3][3] = {0};
	matrix_mul(matrix_1, r_x, r_z);
	matrix_mul(matrix_2, matrix_1, r_y);

	double ret[3];

	matrix_vec_mul(ret, matrix_2, (double[]){src->x, src->y, src->z});

	vec_add(ret, ret, (double[]){
		(N0 * (1 - e_2) + H0) * sin(B0),
		(N0 + H0) * cos(B0) * cos(L0),
		(N0 + H0) * cos(B0) * sin(L0)
	});

	dst->x = ret[0];
	dst->y = ret[1];
	dst->z = ret[2];

	return 0;
}

int lcs_normal_xyz_v(const struct gcs_param *gcs, const struct lcs_param *lcs_param, const struct coord *src_v, struct coord *dst_v)
{
	if (!gcs || !lcs_param || !src_v || !dst_v)
	{
		return -1;
	}
	if (lcs_param->xi || lcs_param->eta)
	{
		return -1;
	}

	double B0 = lcs_param->coord.latitude * ARC_2_DEGREE;
	double L0 = lcs_param->coord.longitude * ARC_2_DEGREE;
	double A = lcs_param->A * ARC_2_DEGREE;

	double r_z[3][3] = {
		{cos(B0), sin(B0), 0},
		{-sin(B0), cos(B0), 0},
		{0, 0, 1}
	};
	double r_x[3][3] = {
		{1, 0, 0},
		{0, cos(-L0), sin(-L0)},
		{0, -sin(-L0), cos(-L0)}
	};
	double r_y[3][3] = {
		{cos(A), 0, sin(A)},
		{0, 1, 0},
		{-sin(A), 0, cos(A)}
	};

	double matrix_1[3][3] = {0}, matrix_2[3][3] = {0};
	matrix_mul(matrix_1, r_x, r_z);
	matrix_mul(matrix_2, matrix_1, r_y);

	double ret[3];

	matrix_vec_mul(ret, matrix_2, (double[]){src_v->dx, src_v->dy, src_v->dz});

	dst_v->dx = ret[0];
	dst_v->dy = ret[1];
	dst_v->dz = ret[2];
	return 0;
}

// 空间球心直角坐标系-> 发射坐标系 (法线坐标系)
int xyz_lcs_normal(const struct gcs_param *gcs_param, const struct lcs_param *lcs_param, const struct coord *src, struct coord *dst)
{
	if (!gcs_param || !lcs_param || !src || !dst)
	{
		return -1;
	}

	if (lcs_param->xi || lcs_param->eta)
	{
		return -2;
	}

	struct coord temp = {0}, lcs_xyz = {0};

	gcs_xyz(gcs_param, &lcs_param->coord, &lcs_xyz);
	temp.x = src->x - lcs_xyz.x;
	temp.y = src->y - lcs_xyz.y;
	temp.z = src->z - lcs_xyz.z;


	double L0 = lcs_param->coord.longitude * ARC_2_DEGREE;
	double B0 = lcs_param->coord.latitude * ARC_2_DEGREE;
	double A = lcs_param->A * ARC_2_DEGREE;

	double matrix_1[3][3] = {0.0}, matrix_2[3][3] = {0.0};

	double r_x[3][3] = {
		{1, 0, 0},
		{0, cos(L0), sin(L0)},
		{0, -sin(L0), cos(L0)}
	};

	double r_y[3][3] = {
		{cos(-A), 0, sin(-A)},
		{0, 1, 0},
		{-sin(-A), 0, cos(-A)}
	};

	double r_z[3][3] = {
		{cos(-B0), sin(-B0), 0},
		{-sin(-B0), cos(-B0), 0},
		{0, 0, 1}
	};

	matrix_mul(matrix_1, r_z, r_x);
	matrix_mul(matrix_2, matrix_1, r_y);
	
	double ret[3] = {0};
	matrix_vec_mul(ret, matrix_2, (double[]){temp.x, temp.y, temp.z});

	dst->x = ret[0];
	dst->y = ret[1];
	dst->z = ret[2];

	return 0;
}

int xyz_lcs_normal_v(const struct gcs_param *gcs_param, const struct lcs_param *lcs_param, const struct coord *src_v, struct coord *dst_v)
{
	if (!gcs_param || !lcs_param || !src_v || !dst_v)
	{
		return -1;
	}
	if (lcs_param->xi || lcs_param->eta)
	{
		return -2;
	}

	double L0 = lcs_param->coord.longitude * ARC_2_DEGREE;
	double B0 = lcs_param->coord.latitude * ARC_2_DEGREE;
	double A = lcs_param->A * ARC_2_DEGREE;

	double matrix_1[3][3] = {0.0}, matrix_2[3][3] = {0.0};

	double r_x[3][3] = {
		{1, 0, 0},
		{0, cos(L0), sin(L0)},
		{0, -sin(L0), cos(L0)}
	};

	double r_y[3][3] = {
		{cos(-A), 0, sin(-A)},
		{0, 1, 0},
		{-sin(-A), 0, cos(-A)}
	};

	double r_z[3][3] = {
		{cos(-B0), sin(-B0), 0},
		{-sin(-B0), cos(-B0), 0},
		{0, 0, 1}
	};

	matrix_mul(matrix_1, r_z, r_x);
	matrix_mul(matrix_2, matrix_1, r_y);
	
	double ret[3] = {0};
	matrix_vec_mul(ret, matrix_2, (double[]){src_v->dx, src_v->dy, src_v->dz});

	dst_v->dx = ret[0];
	dst_v->dy = ret[1];
	dst_v->dz = ret[2];

	return 0;

}

int gcs_coord_transform(const struct coord *src, enum gcs_type src_type, struct coord *dst, enum gcs_type dst_type)
{
	int r;
	struct coord temp;
	switch (src_type)
	{
		case BEI_JING_54:
			r= gcs_xyz(&PARAM_BEI_JING_54, src, &temp);
			break;
		case XI_AN_80:
			r = gcs_xyz(&PARAM_XI_AN_80, src, &temp);
			break;
		case WGS_84:
			r = gcs_xyz(&PARAM_WGS_84, src, &temp);
			break;
		case CGCS_2000:
			r = gcs_xyz(&PARAM_CGCS_2000, src, &temp);
			break;
		default:
			return -2;
	}

	if (r)
		return r;

	switch (dst_type)
	{
		case BEI_JING_54:
			r = xyz_gcs(&PARAM_BEI_JING_54, &temp, dst);
			break;
		case XI_AN_80:
			r = xyz_gcs(&PARAM_XI_AN_80, &temp, dst);
			break;
		case WGS_84:
			r = xyz_gcs(&PARAM_WGS_84, &temp, dst);
			break;
		case CGCS_2000:
			r = xyz_gcs(&PARAM_CGCS_2000, &temp, dst);
			break;
		default:
			return -3;
	}
	return r;
}

// 站心坐标转换（同一参考椭球）
int lcccs_coord_transform(enum gcs_type type, const struct lcccs_param *src_param, const struct coord *src,		\
							const struct lcccs_param *dst_param, struct coord *dst)
{
	return lcccs_coord_transform1(type, src_param, src, type, dst_param, dst);
}

// 站心坐标转换（不同参考椭球）
int lcccs_coord_transform1(enum gcs_type src_type, const struct lcccs_param *src_param, const struct coord *src, \
							enum gcs_type dst_type, const struct lcccs_param *dst_param, struct coord *dst)
{
	int r;
	struct coord temp;
	switch (src_type)
	{
		case BEI_JING_54:
			r = lcccs_normal_xyz(&PARAM_BEI_JING_54, src_param, src, &temp);
			break;
		case XI_AN_80:
			r = lcccs_normal_xyz(&PARAM_XI_AN_80, src_param, src, &temp);
			break;
		case WGS_84:
			r = lcccs_normal_xyz(&PARAM_WGS_84, src_param, src, &temp);
			break;
		case CGCS_2000:
			r = lcccs_normal_xyz(&PARAM_CGCS_2000, src_param, src, &temp);
			break;
		default:
			return -1;
	}
	
	if (r)
		return r;

	switch (dst_type)
	{
		case BEI_JING_54:
			r = xyz_lcccs_normal(&PARAM_BEI_JING_54, dst_param, &temp, dst);
			break;
		case XI_AN_80:
			r = xyz_lcccs_normal(&PARAM_XI_AN_80, dst_param, &temp, dst);
			break;
		case WGS_84:
			r = xyz_lcccs_normal(&PARAM_WGS_84, dst_param, &temp, dst);
			break;
		case CGCS_2000:
			r = xyz_lcccs_normal(&PARAM_CGCS_2000, dst_param, &temp, dst);
			break;
		default:
			return -1;
	}

	return r;
}

int lcccs_lcs_coord_transform(enum gcs_type gcs_type, const struct lcccs_param *src_param, const struct coord *src, const struct lcs_param *dst_param, struct coord *dst)
{
	return lcccs_lcs_coord_transform1(gcs_type, src_param, src, gcs_type, dst_param, dst);
}

int lcccs_lcs_coord_transform1(enum gcs_type src_type, const struct lcccs_param *src_param, const struct coord *src, enum gcs_type dst_type, const struct lcs_param *dst_param, struct coord *dst)
{
	if (!src_param || !src || !dst_param || !dst)
	{
		return -1;
	}
	if (src_param->xi || src_param->eta || dst_param->xi || dst_param->eta)
	{
		return -2;
	}

	struct coord temp = {0};
	int r = 0;

	switch (src_type)
	{
		case BEI_JING_54:
			r = lcccs_normal_xyz(&PARAM_BEI_JING_54, src_param, src, &temp);
			break;
		case XI_AN_80:
			r = lcccs_normal_xyz(&PARAM_XI_AN_80, src_param, src, &temp);
			break;
		case WGS_84:
			r = lcccs_normal_xyz(&PARAM_WGS_84, src_param, src, &temp);
			break;
		case CGCS_2000:
			r = lcccs_normal_xyz(&PARAM_CGCS_2000, src_param, src, &temp);
			break;
		default:
			return -3;
	}
	if (r)
	{
		return r;
	}

	switch (dst_type)
	{
		case BEI_JING_54:
			r = xyz_lcs_normal(&PARAM_BEI_JING_54, dst_param, &temp, dst);
			break;
		case XI_AN_80:
			r = xyz_lcs_normal(&PARAM_XI_AN_80, dst_param, &temp, dst);
			break;
		case WGS_84:
			r = xyz_lcs_normal(&PARAM_WGS_84, dst_param, &temp, dst);
			break;
		case CGCS_2000:
			r = xyz_lcs_normal(&PARAM_CGCS_2000, dst_param, &temp, dst);
			break;
		default:
			return -3;
	}
	return r;
}


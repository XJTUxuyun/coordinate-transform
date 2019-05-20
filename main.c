#include <stdio.h>

#include "coord.h"

const static struct gcs_param PARAM_WGS_84 = {
	.a = 6378137.0,
	.b = 6356752.314245,
	.alpha = 1.0/98.257223563,
};


int main(int argc, char **argv)
{
	printf("a->%f b->%f\n", PARAM_WGS_84.a, PARAM_WGS_84.b);

	printf("----------------test LBH<->XYZ----------------\n");
	struct coord coord1 = {.longitude = 90, .latitude=0, .altitude=10};
	printf("\tlbh_src-> longitude=%f, latitude=%f altitude=%f\n", coord1.longitude, coord1.latitude, coord1.altitude);
	struct coord coord2;
	gcs_xyz(&PARAM_WGS_84, &coord1, &coord2);
	printf("\tlbh->xyz x=%f, y=%f z=%f\n", coord2.longitude, coord2.latitude, coord2.altitude);
	struct coord coord3;
	xyz_gcs(&PARAM_WGS_84, &coord2, &coord3);
	printf("\txyz->lbh longitude=%f, latitude=%f altitude=%f\n", coord3.longitude, coord3.latitude, coord3.altitude);



	printf("----------------test LCCCS<->XYZ----------------\n");
	struct lcccs_param device_1 = {
		.coord = {
			.longitude = 90,
			.latitude = 0,
			.altitude = 0
		},
		.xi = 0.0,
		.eta = 0.0
	};

	struct coord point_1 = {
		.x = 10,
		.y = 20,
		.z = 30
	};

	printf("\tlcccs_src-> x=%f, y=%f z=%f\n", point_1.longitude, point_1.latitude, point_1.altitude);
	struct coord point_1_xyz;

	int r = lcccs_normal_xyz(&PARAM_WGS_84, &device_1, &point_1, &point_1_xyz);

	printf("\tlcccs->xyz status=%d x=%f y=%f z=%f\n", r, point_1_xyz.x, point_1_xyz.y, point_1_xyz.z);

	struct coord point_1_lbh;

	r = xyz_lcccs_normal(&PARAM_WGS_84, &device_1, &point_1_xyz, &point_1_lbh);

	printf("\txyz->lcccs status=%d x=%f y=%f z=%f\n", r, point_1_lbh.longitude, point_1_lbh.latitude, point_1_lbh.altitude);


	printf("----------------test LCS<->XYZ----------------\n");
	struct lcs_param lcs_param = {
		.coord = {
			.longitude = 0,
			.latitude = 0,
			.altitude = 0
		},
		.A = 0.0,
		.xi = 0.0,
		.eta = 0.0
	};

	struct coord lcs_point = {.x = -10, .y = 0, .z = 789};

	printf("\tlcs_src-> x=%f, y=%f z=%f\n", lcs_point.longitude, lcs_point.latitude, lcs_point.altitude);
	struct coord lcs_xyz = {0};
	r = lcs_normal_xyz(&PARAM_WGS_84, &lcs_param, &lcs_point, &lcs_xyz);
	printf("\tlcs->xyz %d x=%f, y=%f, z=%f\n", r, lcs_xyz.x, lcs_xyz.y, lcs_xyz.z);
	struct coord xyz_lcs;
	r = xyz_lcs_normal(&PARAM_WGS_84, &lcs_param, &lcs_xyz, &xyz_lcs);
	printf("\txyz->lcs %d x=%f, y=%f, z=%f\n", r, xyz_lcs.x, xyz_lcs.y, xyz_lcs.z);

	printf("----------------test LCS<->LCCCS----------------\n");
	struct coord lcs_1 = {.x = 123, .y = 456, .z=789};
	printf("\tlcs_src-> x=%f, y=%f z=%f\n", lcs_1.x, lcs_1.y, lcs_1.z);
	struct coord lcs_xyz_1;
	lcs_normal_xyz(&PARAM_WGS_84, &lcs_param, &lcs_1, &lcs_xyz_1);
	struct coord lcccs_xyz;
	xyz_lcccs_normal(&PARAM_WGS_84, &device_1, &lcs_xyz_1, &lcccs_xyz);
	printf("\tlcs->lcccs-> x=%f y=%f z=%f\n", lcccs_xyz.x, lcccs_xyz.y, lcccs_xyz.z);
	struct coord xyz_lcccs;
	lcccs_normal_xyz(&PARAM_WGS_84, &device_1, &lcccs_xyz, &xyz_lcccs);
	struct coord lcs_2;
	xyz_lcs_normal(&PARAM_WGS_84, &lcs_param, &xyz_lcccs, &lcs_2);
	printf("\tlcccs->lcs-> x=%f, y=%f z=%f\n", lcs_2.x, lcs_2.y, lcs_2.z);

	return 0;
}

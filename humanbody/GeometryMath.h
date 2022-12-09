#ifndef _GEOMETRY_MATH_
#define _GEOMETRY_MATH_

#include <math.h>


	namespace Math
	{
		/** 
		* rotate a vertex
		* \param rotation the quaternion that the vertex to rotate 
		* \param vec vertex to rotate
		* \return V result position
		*/
		template <typename V>
		inline V rotate(V & rotation, V & vec) 
		{
			if (rotation == V(0, 0, 0)) 
			{
				return vec;
			}
			double angle = rotation.norm();
			V axis(rotation / angle);
			V res;
			res = vec * cos(angle) + (axis % vec) * sin(angle) + axis * (axis | vec) * (1 - cos(angle));
			return res;

			/*
			const double q0 = cos(angle / 2.0);
			const double sina = sin(angle / 2.0);
			const double q1 = sina * axis[0];
			const double q2 = sina * axis[1];
			const double q3 = sina * axis[2];
			OpenMesh::Vec3d res;
			res[0] = (q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3) * vec[0] + 2 * (q1 * q2 - q0 * q3) * vec[1] + 
			2 * (q1 * q3 + q0 * q2) * vec[2];
			res[1] = 2 * (q2 * q1 + q0 * q3) * vec[0] + (q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3) * vec[1] + 
			2 * (q2 * q3 - q0 * q1) * vec[2];
			res[2] = 2 * (q3 * q1 - q0 * q2) * vec[0] + 2 * (q3 * q2 + q0 * q1) * vec[1] + 
			(q0 * q0 - q1 * q1 - q2 * q2 - q3 * q3) * vec[2];*/
			//Given angle r in radians and unit vector u = ai + bj + ck or [a,b,c]', define:
			//q0 = cos(r/2),  q1 = sin(r/2) a,  q2 = sin(r/2) b,  q3 = sin(r/2) c
			//
			//		[(q0^2 + q1^2 - q2^2 - q3^2)	2(q1q2 - q0q3)				2(q1q3 + q0q2)				]
			//Q  =	[2(q2q1 + q0q3)					(q0^2 - q1^2 + q2^2 - q3^2)	2(q2q3 - q0q1)				]
			//		[2(q3q1 - q0q2)					2(q3q2 + q0q1)				(q0^2 - q1^2 - q2^2 + q3^2)	]
		}

		/** 
		 * judge if a point is on the normal side of a plane
		 * \param v the point on the plane.
		 * \param n the normal of the plane.
		 * \param t the point to be judged
		 * \return bool true or false
		 */
		template <typename V>
		bool normal_side(V v, V n, V t)
		{
			V d;
			d = t - v;
			if ((d | n) < 0)
			{
				return false;
			}
			
			return true;
		}

		/** 
		 * return a vertical vector of n
		 */
		template <typename V>
		V verticalVector(V n)
		{
			/* n[0] * x + n[1] * y + n[2] * z = 0 */
			V v;
			if ((n[0] != 0) && (n[1] == 0) && (n[2] == 0))
			{
				return V(0, 1, 0);
			}
			else if ((n[0] == 0) && (n[1] != 0) && (n[2] == 0))
			{
				return V(0, 0, 1);
			}
			else if ((n[0] == 0) && (n[1] == 0) && (n[2] != 0))
			{
				return V(1, 0, 0);
			}
			else if ((n[0] != 0) && (n[1] != 0) && (n[2] == 0))
			{
				return V(1, -(n[0] / n[1]), 0);
			}
			else if ((n[0] == 0) && (n[1] != 0) && (n[2] != 0))
			{
				return V(0, 1, - (n[1] / n[2]));
			}
			else if ((n[0] != 0) && (n[1] == 0) && (n[2] != 0))
			{
				return V(-(n[2] / n[0]), 0, 1);
			}
			else
			{
				return V(1, 1, -(n[0] + n[1]) / n[2]);
			}
		}

		namespace Geometry
		{
			/**
			* Project an x,y pair onto a sphere of radius r OR a hyperbolic sheet
			* if we are away from the center of the sphere.
			*/
			template <typename T0, typename T1>
			bool project_to_sphere(T0 & v, T1 x, T1 y, T1 width, T1 height)
			{
				if((x >= 0) && (x < width ) &&
				   (y >= 0) && (y < height) )
				{
					double dx;
					double dy;
					double sinx;
					double siny;
					double sinx2siny2;
					dx = (double)(x - 0.5 * width) / (double)width;
					dy = (double)(0.5 * height - y) / (double)height;

					sinx       = sin(M_PI_2 * dx);
					siny       = sin(M_PI_2 * dy);
					sinx2siny2 = sinx * sinx + siny * siny;

					v[0] = sinx;
					v[1] = siny;

					/** 
					* if the point is outside the tracker sphere, set z to 0,
					* else let x, y, z are contented x * x + y * y + z * z = 1
					**/
					v[2] = sinx2siny2 < 1.0 ? sqrt(1 - sinx2siny2) : 0.0;

					return true;
				}
				return false;
			}

			
			/**
			* simulate a track-ball.  Project the points onto the virtual
			* trackball, then figure out the axis of rotation, which is the cross
			* product of v0 v1 and O v0 (O is the center of the ball, 0,0,0)
			* Note:  This is a deformed trackball-- is a trackball in the center,
			* but is deformed into a hyperbolic sheet of rotation away from the
			* center.  This particular function was chosen after trying out
			* several variations.
			*
			* It is assumed that the arguments to this routine are in the range
			* (-1.0 ... 1.0)
			*/
			template <typename T0, typename T1, typename T2>
			void rotation_vector(T0 & v, T1 & x0, T1 & y0, T1 &  x1, T1 & y1, T2 & width, T2 & height)
			{
				v.zero();
				T0 v0;
				T0 v1;
				if (project_to_sphere(v0, x0, y0, width, height) && project_to_sphere(v1, x1, y1, width, height))
				{
					v = v0 % v1;
					v.normalize();
					double cos_angle;
					cos_angle = (v0 | v1);
					double angle;
					angle = acos(std::min(1.0, std::max(cos_angle, -1.0))); 
					v *= angle;
				}//2.0 *  * 180.0 / M_PI
			}


			template<typename MatrixType, typename VectorType> 
		MatrixType rotationMatrix(const VectorType& unitVec, double angle)
		{
			//angle = -angle;
			
			double sin_a = std::sin( angle / 2 );
			double cos_a = std::cos( angle / 2 );
			double X = unitVec[0] * sin_a;
			double Y = unitVec[1] * sin_a;
			double Z = unitVec[2] * sin_a;
			double W = cos_a;

			double xx = X * X;
			double xy  = X * Y;
			double xz  = X * Z;
			double xw = X * W;

			double yy = Y * Y;
			double yz = Y * Z;
			double yw = Y * W;

			double zz = Z * Z;
			double zw = Z * W;

			MatrixType resMatrix;
			/*resMatrix << 1. - 2. * ( yy + zz ), 2. * ( xy + zw ),      2. * ( xz - yw ),    
									2. * ( xy - zw ),      1. - 2. * ( xx + zz ), 2. * ( yz + xw ), 
									2. * ( xz + yw ),      2. * ( yz - xw ),      1. - 2. * ( xx + yy );*/
			resMatrix.col(0) << 1. - 2. * ( yy + zz ), 2. * ( xy - zw ),  2. * ( xz + yw );
			resMatrix.col(1) << 2. * ( xy + zw ), 1. - 2. * ( xx + zz ),  2. * ( yz - xw ); 
			resMatrix.col(2) << 2. * ( xz - yw ), 2. * ( yz + xw ), 1. - 2. * ( xx + yy );
			
			return resMatrix;
		}
		
		template<typename MatrixType, typename VectorType> 
		MatrixType rotationMatrix(const VectorType& Vec)
		{
			double angle = Vec.norm();
			VectorType vec;
			if (abs(angle) > 0.00001)
			{
				vec = Vec / angle;
			}
			else
			{
				vec = VectorType(0., 0., 0.);
				angle = 0.;
			}
			return rotationMatrix<MatrixType, VectorType>(vec, angle);
		}

		template<typename MatrixType, typename VectorType>
		void Matrix2Quat(const MatrixType m, VectorType& q)
		{
			#define QX q[0]
			#define QY q[1]
			#define QZ q[2]
			#define QW q[3]
			
			//const MatrixType m = _m.transpose();

			double s;
			double tq[4];
			int    i, j;

			// Use tq to store the largest trace
			tq[0] = 1 + m(0, 0)+m(1, 1)+m(2, 2);
			tq[1] = 1 + m(0, 0)-m(1, 1)-m(2, 2);
			tq[2] = 1 - m(0, 0)+m(1, 1)-m(2, 2);
			tq[3] = 1 - m(0, 0)-m(1, 1)+m(2, 2);

			// Find the maximum (could also use stacked if's later)
			j = 0;
			for(i=1;i<4;i++) j = (tq[i]>tq[j])? i : j;

			// check the diagonal
			if (j==0)
			{
				/* perform instant calculation */
				QW = tq[0];
				QX = m(1, 2)-m(2, 1);
				QY = m(2, 0)-m(0, 2);
				QZ = m(0, 1)-m(1, 0);
			}
			else if (j==1)
			{
				QW = m(1, 2)-m(2, 1);
				QX = tq[1];
				QY = m(0, 1)+m(1, 0);
				QZ = m(2, 0)+m(0, 2);
			}
			else if (j==2)
			{
				QW = m(2, 0)-m(0, 2);
				QX = m(0, 1)+m(1, 0);
				QY = tq[2];
				QZ = m(1, 2)+m(2, 1);
			}
			else /* if (j==3) */
			{
				QW = m(0, 1)-m(1, 0);
				QX = m(2, 0)+m(0, 2);
				QY = m(1, 2)+m(2, 1);
				QZ = tq[3];
			}

			s = sqrt(0.25/tq[j]);
			QW *= s;
			QX *= s;
			QY *= s;
			QZ *= s;
		}
		template<typename VectorType0, typename VectorType1>
		void axisAndAngle(const VectorType0& q, VectorType1& rot_axis)
		{
			double sin_theta_half = sqrt(1 - q[3] * q[3]);
			
			rot_axis[0] = q[0];
			rot_axis[1] = q[1];
			rot_axis[2] = q[2];
			
			if (abs(sin_theta_half) > 0.00001)
			{
				rot_axis[0] /= sin_theta_half;
				rot_axis[1] /= sin_theta_half;
				rot_axis[2] /= sin_theta_half;
			}
			else
			{
				sin_theta_half = 0.;
			}
			
			double rot_angle = asin(sin_theta_half) * 2;

			if (q[3] < 0 && rot_angle != 0.)
			{
				/*rot_angle = 2 * (M_PI - rot_angle);
				rot_angle = 2 * M_PI - rot_angle;*/
				//rot_angle *= 2; 
				rot_axis = -rot_axis;
			}
			
			rot_axis *= rot_angle;
		}
		}

	}

#endif
/*
Copyright 2011 Frans van den Bergh. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of
      conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list
      of conditions and the following disclaimer in the documentation and/or other materials
      provided with the distribution.

THIS SOFTWARE IS PROVIDED BY Frans van den Bergh ''AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL Frans van den Bergh OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those of the
authors and should not be interpreted as representing official policies, either expressed
or implied, of the Council for Scientific and Industrial Research (CSIR).
*/

#ifndef ECCENTRICITY_CORRECTOR_H
#define ECCENTRICITY_CORRECTOR_H

// An implementation of
// Ahn, S.J., Warnecke, H.-J., Systematic geometric image measurement errors of circular 
// target objects: mathematical formulation and correction, Photogrammetric Record,
// 16(93), 485-502, 1999

#include <Eigen/Dense>

class Eccentricity_corrector {
  public:
    Eccentricity_corrector(const Eigen::Matrix3d& rotation, const Eigen::Vector3d& cop,
        double circle_radius, double focal_length) 
    : R_star(rotation), N_p(0,0,1), X_star(cop), r(circle_radius), c(focal_length)
    {
        x = N_p.cross(R_star.row(2));
        
        x /= x.norm(); // TODO: check for div-by-zero
        
        z = N_p;
        y = z.cross(x);
        y /= y.norm(); // TODO: check for divide-by-zero
        
        R_o.row(0) = x;
        R_o.row(1) = y;
        R_o.row(2) = z;
        
        oR_star = R_star * (R_o.transpose());
        
        // under the assumption that phi (the yaw angle) is zero, 
        // we know that elements (2,2) and (0,0) of oR_star are formed
        // from only a single cosine (i.e., the phi rotation matrix has no effect on them)
        
        // omega is the tilt of the camera relative to target plane
        omega = acos(oR_star(2,2));
        
        // kappa is azimuth angle between camera and target plane 
        kappa = acos(oR_star(0,0));
    }
    
    // this eccentricity vector must be subtracted from the normalized image
    // coordinates of the ellipse centre
    Eigen::Vector2d eccentricity(const Eigen::Vector3d& world_point) const {
        Eigen::Vector3d X_p = world_point;
        
        double delta_y = y.dot(X_star - X_p);
        double delta_z = z.dot(X_star - X_p);
        double alpha = atan2(-delta_y, delta_z);
        double d = sqrt(SQR(delta_y) + SQR(delta_z));
        double x_p = (X_p - X_star).dot(x);
        
        double l = d * cos(omega - alpha);
        double so = sin(omega);
        double eup = (so*so*c*x_p/l) / (SQR(l/r) - so*so);
        double evp = (-c*(d/l)*so*cos(alpha))/(SQR(l/r) - so*so);
        
        double ck = cos(kappa);
        double sk = sin(kappa);
        double eu = ck*eup + sk*evp;
        double ev = -(-sk*eup + ck*evp); // flip y-axis, because we use top left (0,0) convention
        
        return Eigen::Vector2d(eu, ev);
    }
    
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    Eigen::Matrix3d R_star;        // camera rotation matrix
    Eigen::Vector3d N_p;           // unit normal vector of target in world coords
    Eigen::Vector3d x;             // basis vector of target coordinate system
    Eigen::Vector3d y;             // basis vector of target coordinate system
    Eigen::Vector3d z;             // basis vector of target coordinate system
    Eigen::Matrix3d R_o;           
    Eigen::Matrix3d oR_star;       // rotation matrix going from intermediate target frame to intermediate camera frame
    Eigen::Vector3d X_star;        // COP in world frame
    
    double omega;
    double kappa;
    double r;
    double c;
};

#endif

#!/usr/bin/env python

# Copyright (C) 2017 Electric Movement Inc.
#
# This file is part of Robotic Arm: Pick and Place project for Udacity
# Robotics nano-degree program
#
# All Rights Reserved.

# Author: Harsh Pandya

# import modules
import rospy
import tf
from kuka_arm.srv import *
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
from geometry_msgs.msg import Pose
from mpmath import *
from sympy import *


def handle_calculate_IK(req):
    rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
    if len(req.poses) < 1:
        print "No valid poses received"
        return -1
    else:
        # Create symbols
        d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8') #link offset
        a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7') #link length
        alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7') # Twist angle
        q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8') # theta_i
        r, p, y = symbols('r p y') # roll, pitch, and yaw symbols

        # Create Modified DH parameters
        dh = {alpha0:     0,   a0:      0,   d1:  0.75,
              alpha1: -pi/2,   a1:   0.35,   d2:     0,   q2: q2-pi/2,
              alpha2:     0,   a2:   1.25,   d3:     0,
              alpha3: -pi/2,   a3: -0.054,   d4:   1.5,
              alpha4:  pi/2,   a4:      0,   d5:     0,
              alpha5: -pi/2,   a5:      0,   d6:     0,
              alpha6:     0,   a6:      0,   d7: 0.303,   q7:       0}

    	#
    	# Define Modified DH Transformation matrix
    	#
    	def TF_Matrix(q, alpha, a, d):
        	T = Matrix([[           cos(q),           -sin(q),           0,             a],
                	    [sin(q)*cos(alpha), cos(q)*cos(alpha), -sin(alpha), -sin(alpha)*d],
                   	    [sin(q)*sin(alpha), cos(q)*sin(alpha),  cos(alpha),  cos(alpha)*d],
                    	    [                0,                 0,           0,             1]])
        	return(T)       

        # Create individual transformation matrices
        T0_1 = TF_Matrix(q1, alpha0, a0, d1).subs(dh)
        T1_2 = TF_Matrix(q2, alpha1, a1, d2).subs(dh)
        T2_3 = TF_Matrix(q3, alpha2, a2, d3).subs(dh)
        T3_4 = TF_Matrix(q4, alpha3, a3, d4).subs(dh)
        T4_5 = TF_Matrix(q5, alpha4, a4, d5).subs(dh)
        T5_6 = TF_Matrix(q6, alpha5, a5, d6).subs(dh)
        T6_EE = TF_Matrix(q7, alpha6, a6, d7).subs(dh)

        # Total Matrix
        T0_EE = T0_1 * T1_2 * T2_3 * T3_4 * T4_5 * T5_6 * T6_EE

    	# Find EE rotation Matrix
    	# Define RPY rotation matrices

    	r, p, y = symbols('r p y')

    	Rot_x = Matrix([[        1,         0,         0],
        	        [        0,    cos(r),   -sin(r)],
                	[        0,    sin(r),    cos(r)]]) # ROLL

    	Rot_y = Matrix([[   cos(p),         0,    sin(p)],
                        [        0,         1,         0],
                        [  -sin(p),         0,    cos(p)]]) # PITCH

    	Rot_z = Matrix([[   cos(y),   -sin(y),         0],
                        [   sin(y),    cos(y),         0],
                        [        0,         0,         1]]) # ROLL

    	Rot_EE = Rot_z * Rot_y * Rot_x

        
	# Extract rotation matrices from the transformation matrices
        # R_corr = R_z(pi) * R_y(-pi/2)
        # R_rpy = R_z(y) * R_y(p) * R_x(r)

        # Initialize service response
        joint_trajectory_list = []
        for x in xrange(0, len(req.poses)):
            # IK code starts here
            joint_trajectory_point = JointTrajectoryPoint()

	    # Extract end-effector position and orientation from request
	    # px,py,pz = end-effector position
	    # roll, pitch, yaw = end-effector orientation
            px = req.poses[x].position.x
            py = req.poses[x].position.y
            pz = req.poses[x].position.z

            (roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
                [req.poses[x].orientation.x, req.poses[x].orientation.y,
                    req.poses[x].orientation.z, req.poses[x].orientation.w])

    	    # More information found in KR210 Forward Kinematics section
    	    Rot_Error = Rot_z.subs(y, radians(180)) * Rot_y.subs(p, radians(-90))

    	    Rot_EE = Rot_EE * Rot_Error
    	    Rot_EE = Rot_EE.subs({'r': roll, 'p': pitch, 'y': yaw})

    	    EE = Matrix([[px],
                         [py],
                         [pz]])

    	    WC = EE - (0.303) * Rot_EE[:, 2]

    	    # Calculate joint angle by using Geometric IK method
    	    # more information found in the inverse Kinematics with Kuka KR210
    	    theta1 = atan2(WC[1], WC[0])

    	    # SSS triangle for theta2 and theta3
    	    side_a = 1.501
    	    side_b = sqrt(pow((sqrt(WC[0]*WC[0] + WC[1] * WC[1]) - 0.35), 2) + pow((WC[2] - 0.75), 2))
    	    side_c = 1.25

    	    angle_a = acos((side_b * side_b + side_c * side_c - side_a * side_a) / (2 * side_b * side_c))
    	    angle_b = acos((side_a * side_a + side_c * side_c - side_b * side_b) / (2 * side_a * side_c))
    	    angle_c = acos((side_a * side_a + side_b * side_b - side_c * side_c) / (2 * side_a * side_b))

    	    theta2 = pi/2 - angle_a -atan2(WC[2] - 0.75, sqrt(WC[0] * WC[0] + WC[1] * WC[1]) -0.35) 
    	    theta3 = pi/2 - (angle_b + 0.036)   # 0.036 accounts for sag in link 4 of  -0.054m

    	    R0_3 = T0_1[0:3, 0:3] * T1_2[0:3, 0:3] * T2_3[0:3, 0:3]
    	    R0_3 = R0_3.evalf(subs={q1: theta1, q2: theta2, q3: theta3})

    	    R3_6 = R0_3.inv("LU") * Rot_EE

    	    # Euler angles from rotation matrix
    	    theta4 = atan2(R3_6[2,2], -R3_6[0,2])
    	    theta5 = atan2(sqrt(R3_6[0,2]*R3_6[0,2] + R3_6[2,2]*R3_6[2,2]), R3_6[1,2])
    	    theta6 = atan2(-R3_6[1,1], R3_6[1,0])



            # Populate response for the IK request
            # In the next line replace theta1,theta2...,theta6 by your joint angle variables
            joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5, theta6]
            joint_trajectory_list.append(joint_trajectory_point)

        rospy.loginfo("length of Joint Trajectory List: %s" % len(joint_trajectory_list))
        return CalculateIKResponse(joint_trajectory_list)


def IK_server():
    # initialize node and declare calculate_ik service
    rospy.init_node('IK_server')
    s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
    print "Ready to receive an IK request"
    rospy.spin()

if __name__ == "__main__":
    IK_server()

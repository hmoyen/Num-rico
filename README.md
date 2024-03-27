# Inertial Measurement Unit (IMU) Data Processing

This repository contains code and documentation for processing data from an Inertial Measurement Unit (IMU) sensor to estimate the state of a moving object. The project utilizes numerical techniques such as numerical integration, Runge-Kutta method, numerical quadrature, and cubic spline interpolation to estimate the position, velocity, and orientation of the object based on IMU data.

## Table of Contents

- [Introduction](#introduction)
- [Numerical Methodology](#numerical-methodology)
  - [Runge-Kutta](#runge-kutta)
  - [Numerical Quadrature](#numerical-quadrature)
  - [Cubic Spline Interpolation](#cubic-spline-interpolation)

## Introduction

The project aims to estimate the state of an object in motion using data collected from an Inertial Measurement Unit (IMU). IMU data typically includes measurements of linear acceleration, angular velocity, and sometimes magnetic field orientation. By processing this data using numerical techniques, we can estimate various parameters such as position, velocity, and orientation of the object.

## Numerical Methodology

To estimate the state of the object, we employ three main numerical techniques:

### Runge-Kutta

The Runge-Kutta method is used to numerically integrate the equations of motion based on the measured accelerations. This method provides a weighted average of linear accelerations, which is then used to update the position and velocity of the object over time.

### Numerical Quadrature

Numerical quadrature is applied to approximate the integral of angular velocity over time. This integral provides the angular displacement, which is crucial for estimating the orientation of the object. We use methods such as the Trapezoidal Rule to approximate this integral.

### Cubic Spline Interpolation

Cubic spline interpolation is utilized to approximate the position function based on the estimated positions obtained from the Runge-Kutta method. This interpolation technique provides a continuous representation of the object's trajectory, allowing for smoother visualization and analysis.


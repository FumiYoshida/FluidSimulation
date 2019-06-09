using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Alea;
using Alea.Parallel;

namespace FluidSimulation
{
    class Mps
    {
        /*=====================================================================
  mps.c   
  (c) Kazuya SHIBATA, Kohei MUROTANI and Seiichi KOSHIZUKA (2014) 

   Fluid Simulation Program Based on a Particle Method (the MPS method)
   Last update: May 21, 2014
=======================================================================*/


        int DIM;
        double PARTICLE_DISTANCE;
        double DT;
        int OUTPUT_INTERVAL;


        int ARRAY_SIZE;
        double FINISH_TIME;
        double KINEMATIC_VISCOSITY;
        double FLUID_DENSITY;
        double G_X;
        double G_Y;
        double G_Z;
        double RADIUS_FOR_NUMBER_DENSITY;
        double RADIUS_FOR_GRADIENT;
        double RADIUS_FOR_LAPLACIAN;
        double COLLISION_DISTANCE;
        double THRESHOLD_RATIO_OF_NUMBER_DENSITY;
        double COEFFICIENT_OF_RESTITUTION;
        double COMPRESSIBILITY;
        double EPS;
        int ON;
        int OFF;
        double RELAXATION_COEFFICIENT_FOR_PRESSURE;
        int GHOST;
        int FLUID;
        int WALL;
        int DUMMY_WALL;
        int GHOST_OR_DUMMY;
        int SURFACE_PARTICLE;
        int INNER_PARTICLE;
        int DIRICHLET_BOUNDARY_IS_NOT_CONNECTED;
        int DIRICHLET_BOUNDARY_IS_CONNECTED;
        int DIRICHLET_BOUNDARY_IS_CHECKED;

        static double[] Acceleration;
        static int[] ParticleType;
        static double[] Position;
        static double[] Velocity;
        static double[] Pressure;
        static double[] NumberDensity;
        static int[] BoundaryCondition;
        static double[] SourceTerm;
        static int[] FlagForCheckingBoundaryCondition;
        static double[] CoefficientMatrix;
        static double[] MinimumPressure;
        double Time;
        int NumberOfParticles;
        double Re_forNumberDensity, Re2_forNumberDensity;
        double Re_forGradient, Re2_forGradient;
        double Re_forLaplacian, Re2_forLaplacian;
        double N0_forNumberDensity;
        double N0_forGradient;
        double N0_forLaplacian;
        double Lambda;
        double collisionDistance, collisionDistance2;
        double FluidDensity;

        List<double[]> positions;
        List<int[]> particletypes;
        int[,] ijindex;
        System.Windows.Forms.DataVisualization.Charting.Chart chart;
        Action DoEvent;
        Gpu gpu;
        bool usegpu;
        public System.Windows.Forms.Label debuglabel;

        public Mps(System.Windows.Forms.DataVisualization.Charting.Chart chart, Action DoEvent, bool usegpu)
        {
            DIM = 2;
            PARTICLE_DISTANCE = 0.025;
            DT = 0.003;
            OUTPUT_INTERVAL = 10;

            /* for three-dimensional simulation */
            /*
            #define DIM                  3
            #define PARTICLE_DISTANCE    0.075
            #define DT                   0.003
            #define OUTPUT_INTERVAL      2 
            */

            ARRAY_SIZE = 5000;
            FINISH_TIME = 20.0;
            KINEMATIC_VISCOSITY = (1.0E-6);
            FLUID_DENSITY = 1000.0;
            G_X = 0;
            G_Y = -9.8;
            G_Z = 0;
            RADIUS_FOR_NUMBER_DENSITY = (2.1 * PARTICLE_DISTANCE);
            RADIUS_FOR_GRADIENT = (2.1 * PARTICLE_DISTANCE);
            RADIUS_FOR_LAPLACIAN = (3.1 * PARTICLE_DISTANCE);
            COLLISION_DISTANCE = (0.5 * PARTICLE_DISTANCE);
            THRESHOLD_RATIO_OF_NUMBER_DENSITY = 0.97;
            COEFFICIENT_OF_RESTITUTION = 0.2;
            COMPRESSIBILITY = (0.45E-9);
            EPS = (0.01 * PARTICLE_DISTANCE);
            ON = 1;
            OFF = 0;
            RELAXATION_COEFFICIENT_FOR_PRESSURE = 0.2;
            GHOST = -1;
            FLUID = 0;
            WALL = 2;
            DUMMY_WALL = 3;
            GHOST_OR_DUMMY = -1;
            SURFACE_PARTICLE = 1;
            INNER_PARTICLE = 0;
            DIRICHLET_BOUNDARY_IS_NOT_CONNECTED = 0;
            DIRICHLET_BOUNDARY_IS_CONNECTED = 1;
            DIRICHLET_BOUNDARY_IS_CHECKED = 2;

            Acceleration = new double[3 * ARRAY_SIZE];
            ParticleType = new int[ARRAY_SIZE];
            Position = new double[3 * ARRAY_SIZE];
            Velocity = new double[3 * ARRAY_SIZE];
            Pressure = new double[ARRAY_SIZE];
            NumberDensity = new double[ARRAY_SIZE];
            BoundaryCondition = new int[ARRAY_SIZE];
            SourceTerm = new double[ARRAY_SIZE];
            FlagForCheckingBoundaryCondition = new int[ARRAY_SIZE];
            CoefficientMatrix = new double[ARRAY_SIZE * ARRAY_SIZE];
            MinimumPressure = new double[ARRAY_SIZE];

            positions = new List<double[]>();
            particletypes = new List<int[]>();
            this.chart = chart;
            this.DoEvent = DoEvent;
            gpu = Gpu.Default;
            this.usegpu = usegpu;

        }


        public void Main()
        {

            Console.WriteLine("\n*** START PARTICLE-SIMULATION ***\n");
            if (DIM == 2)
            {
                initializeParticlePositionAndVelocity_for2dim();
            }
            else
            {
                initializeParticlePositionAndVelocity_for3dim();
            }
            calConstantParameter();
            mainLoopOfSimulation();
            Console.WriteLine("*** END ***\n\n");
        }


        void initializeParticlePositionAndVelocity_for2dim()
        {
            int iX, iY;
            int nX, nY;
            double x, y, z;
            int i = 0;
            int flagOfParticleGeneration;

            nX = (int)(1.0 / PARTICLE_DISTANCE) + 5;
            nY = (int)(0.6 / PARTICLE_DISTANCE) + 5;
            for (iX = -4; iX < nX; iX++)
            {
                for (iY = -4; iY < nY; iY++)
                {
                    x = PARTICLE_DISTANCE * (double)(iX);
                    y = PARTICLE_DISTANCE * (double)(iY);
                    z = 0.0;
                    flagOfParticleGeneration = OFF;

                    /* dummy wall region */
                    if (((x > -4.0 * PARTICLE_DISTANCE + EPS) && (x <= 1.00 + 4.0 * PARTICLE_DISTANCE + EPS)) && ((y > 0.0 - 4.0 * PARTICLE_DISTANCE + EPS) && (y <= 0.6 + EPS)))
                    {
                        ParticleType[i] = DUMMY_WALL;
                        flagOfParticleGeneration = ON;
                    }

                    /* wall region */
                    if (((x > -2.0 * PARTICLE_DISTANCE + EPS) && (x <= 1.00 + 2.0 * PARTICLE_DISTANCE + EPS)) && ((y > 0.0 - 2.0 * PARTICLE_DISTANCE + EPS) && (y <= 0.6 + EPS)))
                    {
                        ParticleType[i] = WALL;
                        flagOfParticleGeneration = ON;
                    }

                    /* wall region */
                    if (((x > -4.0 * PARTICLE_DISTANCE + EPS) && (x <= 1.00 + 4.0 * PARTICLE_DISTANCE + EPS)) && ((y > 0.6 - 2.0 * PARTICLE_DISTANCE + EPS) && (y <= 0.6 + EPS)))
                    {
                        ParticleType[i] = WALL;
                        flagOfParticleGeneration = ON;
                    }

                    /* empty region */
                    if (((x > 0.0 + EPS) && (x <= 1.00 + EPS)) && (y > 0.0 + EPS))
                    {
                        flagOfParticleGeneration = OFF;
                    }

                    /* fluid region */
                    if (((x > 0.0 + EPS) && (x <= 0.25 + EPS)) && ((y > 0.0 + EPS) && (y <= 0.50 + EPS)))
                    {
                        ParticleType[i] = FLUID;
                        flagOfParticleGeneration = ON;
                    }

                    if (flagOfParticleGeneration == ON)
                    {
                        Position[i * 3] = x; Position[i * 3 + 1] = y; Position[i * 3 + 2] = z;
                        i++;
                    }
                }
            }
            NumberOfParticles = i;
            for (i = 0; i < NumberOfParticles * 3; i++) { Velocity[i] = 0.0; }

            ijindex = new int[NumberOfParticles * (NumberOfParticles - 1) / 2, 2];
            int serial = 0;
            for (i = 0; i < NumberOfParticles - 1; i++)
            {
                for (int j = i + 1; j < NumberOfParticles; j++)
                {
                    ijindex[serial, 0] = i;
                    ijindex[serial, 1] = j;
                    serial++;
                }
            }
        }


        void initializeParticlePositionAndVelocity_for3dim()
        {
            int iX, iY, iZ;
            int nX, nY, nZ;
            double x, y, z;
            int i = 0;
            int flagOfParticleGeneration;

            nX = (int)(1.0 / PARTICLE_DISTANCE) + 5;
            nY = (int)(0.6 / PARTICLE_DISTANCE) + 5;
            nZ = (int)(0.3 / PARTICLE_DISTANCE) + 5;
            for (iX = -4; iX < nX; iX++)
            {
                for (iY = -4; iY < nY; iY++)
                {
                    for (iZ = -4; iZ < nZ; iZ++)
                    {
                        x = PARTICLE_DISTANCE * iX;
                        y = PARTICLE_DISTANCE * iY;
                        z = PARTICLE_DISTANCE * iZ;
                        flagOfParticleGeneration = OFF;

                        /* dummy wall region */
                        if ((((x > -4.0 * PARTICLE_DISTANCE + EPS) && (x <= 1.00 + 4.0 * PARTICLE_DISTANCE + EPS)) && ((y > 0.0 - 4.0 * PARTICLE_DISTANCE + EPS) && (y <= 0.6 + EPS))) && ((z > 0.0 - 4.0 * PARTICLE_DISTANCE + EPS) && (z <= 0.3 + 4.0 * PARTICLE_DISTANCE + EPS)))
                        {
                            ParticleType[i] = DUMMY_WALL;
                            flagOfParticleGeneration = ON;
                        }

                        /* wall region */
                        if ((((x > -2.0 * PARTICLE_DISTANCE + EPS) && (x <= 1.00 + 2.0 * PARTICLE_DISTANCE + EPS)) && ((y > 0.0 - 2.0 * PARTICLE_DISTANCE + EPS) && (y <= 0.6 + EPS))) && ((z > 0.0 - 2.0 * PARTICLE_DISTANCE + EPS) && (z <= 0.3 + 2.0 * PARTICLE_DISTANCE + EPS)))
                        {
                            ParticleType[i] = WALL;
                            flagOfParticleGeneration = ON;
                        }

                        /* wall region */
                        if ((((x > -4.0 * PARTICLE_DISTANCE + EPS) && (x <= 1.00 + 4.0 * PARTICLE_DISTANCE + EPS)) && ((y > 0.6 - 2.0 * PARTICLE_DISTANCE + EPS) && (y <= 0.6 + EPS))) && ((z > 0.0 - 4.0 * PARTICLE_DISTANCE + EPS) && (z <= 0.3 + 4.0 * PARTICLE_DISTANCE + EPS)))
                        {
                            ParticleType[i] = WALL;
                            flagOfParticleGeneration = ON;
                        }

                        /* empty region */
                        if ((((x > 0.0 + EPS) && (x <= 1.00 + EPS)) && (y > 0.0 + EPS)) && ((z > 0.0 + EPS) && (z <= 0.3 + EPS)))
                        {
                            flagOfParticleGeneration = OFF;
                        }

                        /* fluid region */
                        if ((((x > 0.0 + EPS) && (x <= 0.25 + EPS)) && ((y > 0.0 + EPS) && (y < 0.5 + EPS))) && ((z > 0.0 + EPS) && (z <= 0.3 + EPS)))
                        {
                            ParticleType[i] = FLUID;
                            flagOfParticleGeneration = ON;
                        }

                        if (flagOfParticleGeneration == ON)
                        {
                            Position[i * 3] = x;
                            Position[i * 3 + 1] = y;
                            Position[i * 3 + 2] = z;
                            i++;
                        }
                    }
                }
            }
            NumberOfParticles = i;
            for (i = 0; i < NumberOfParticles * 3; i++) { Velocity[i] = 0.0; }
        }


        void calConstantParameter()
        {

            Re_forNumberDensity = RADIUS_FOR_NUMBER_DENSITY;
            Re_forGradient = RADIUS_FOR_GRADIENT;
            Re_forLaplacian = RADIUS_FOR_LAPLACIAN;
            Re2_forNumberDensity = Re_forNumberDensity * Re_forNumberDensity;
            Re2_forGradient = Re_forGradient * Re_forGradient;
            Re2_forLaplacian = Re_forLaplacian * Re_forLaplacian;
            calNZeroAndLambda();
            FluidDensity = FLUID_DENSITY;
            collisionDistance = COLLISION_DISTANCE;
            collisionDistance2 = collisionDistance * collisionDistance;
            Time = 0.0;
        }


        void calNZeroAndLambda()
        {
            int iZ_start, iZ_end;
            if (DIM == 2)
            {
                iZ_start = 0; iZ_end = 1;
            }
            else
            {
                iZ_start = -4; iZ_end = 5;
            }

            N0_forNumberDensity = 0.0;
            N0_forGradient = 0.0;
            N0_forLaplacian = 0.0;
            Lambda = 0.0;
            double xi = 0.0;
            double yi = 0.0;
            double zi = 0.0;

            for (int iX = -4; iX < 5; iX++)
            {
                for (int iY = -4; iY < 5; iY++)
                {
                    for (int iZ = iZ_start; iZ < iZ_end; iZ++)
                    {
                        if (((iX == 0) && (iY == 0)) && (iZ == 0)) continue;
                        double xj = PARTICLE_DISTANCE * (double)(iX);
                        double yj = PARTICLE_DISTANCE * (double)(iY);
                        double zj = PARTICLE_DISTANCE * (double)(iZ);
                        double distance2 = (xj - xi) * (xj - xi) + (yj - yi) * (yj - yi) + (zj - zi) * (zj - zi);
                        double distance = Math.Sqrt(distance2);
                        N0_forNumberDensity += weight(distance, Re_forNumberDensity);
                        N0_forGradient += weight(distance, Re_forGradient);
                        N0_forLaplacian += weight(distance, Re_forLaplacian);
                        Lambda += distance2 * weight(distance, Re_forLaplacian);
                    }
                }
            }
            Lambda = Lambda / N0_forLaplacian;
        }


        double weight(double distance, double re)
        {
            double weightIJ;

            if (distance >= re)
            {
                weightIJ = 0.0;
            }
            else
            {
                weightIJ = (re / distance) - 1.0;
            }
            return weightIJ;
        }


        void mainLoopOfSimulation()
        {
            int iTimeStep = 0;

            // WriteData();
            ShowData();

            while (true)
            {
                calGravity();
                calViscosity();
                moveParticle();
                collision();
                calPressure();
                calPressureGradient();
                moveParticleUsingPressureGradient();

                iTimeStep++;
                Time += DT;
                debuglabel.Text = "time: " + Time + "second";
                if ((iTimeStep % OUTPUT_INTERVAL) == 0)
                {
                    Console.WriteLine("TimeStepNumber: %4d   Time: %lf(s)   NumberOfParticless: %d\n", iTimeStep, Time, NumberOfParticles);
                    // WriteData();
                    ShowData();
                }
                if (Time >= FINISH_TIME) { break; }
            }
        }


        void calGravity()
        {
            int i;

            for (i = 0; i < NumberOfParticles; i++)
            {
                if (ParticleType[i] == FLUID)
                {
                    Acceleration[i * 3] = G_X;
                    Acceleration[i * 3 + 1] = G_Y;
                    Acceleration[i * 3 + 2] = G_Z;
                }
                else
                {
                    Acceleration[i * 3] = 0.0;
                    Acceleration[i * 3 + 1] = 0.0;
                    Acceleration[i * 3 + 2] = 0.0;
                }
            }
        }


        void calViscosity()
        {
            double a = (KINEMATIC_VISCOSITY) * (2.0 * DIM) / (N0_forLaplacian * Lambda);
            if (!usegpu)
            {
                for (int i = 0; i < NumberOfParticles; i++)
                {
                    if (ParticleType[i] != FLUID) continue;
                    double viscosityTerm_x = 0.0;
                    double viscosityTerm_y = 0.0;
                    double viscosityTerm_z = 0.0;

                    for (int j = 0; j < NumberOfParticles; j++)
                    {
                        if ((j == i) || (ParticleType[j] == GHOST)) continue;
                        double xij = Position[j * 3] - Position[i * 3];
                        double yij = Position[j * 3 + 1] - Position[i * 3 + 1];
                        double zij = Position[j * 3 + 2] - Position[i * 3 + 2];
                        double distance2 = (xij * xij) + (yij * yij) + (zij * zij);
                        double distance = Math.Sqrt(distance2);
                        if (distance < Re_forLaplacian)
                        {
                            double w = weight(distance, Re_forLaplacian);
                            viscosityTerm_x += (Velocity[j * 3] - Velocity[i * 3]) * w;
                            viscosityTerm_y += (Velocity[j * 3 + 1] - Velocity[i * 3 + 1]) * w;
                            viscosityTerm_z += (Velocity[j * 3 + 2] - Velocity[i * 3 + 2]) * w;
                        }
                    }
                    viscosityTerm_x = viscosityTerm_x * a;
                    viscosityTerm_y = viscosityTerm_y * a;
                    viscosityTerm_z = viscosityTerm_z * a;
                    Acceleration[i * 3] += viscosityTerm_x;
                    Acceleration[i * 3 + 1] += viscosityTerm_y;
                    Acceleration[i * 3 + 2] += viscosityTerm_z;
                }
            }
            else
            {
                // GPUを使うとき
                // Aleaではクラスの変数が使えないので、関数内の変数にしてから使う
                // また、weightの関数を呼び出さず、直接四則演算の形で書く。
                var particletype = ParticleType;
                var position = Position;
                var acceleration = Acceleration;
                var velocity = Velocity;
                var re = Re2_forLaplacian;
                var fluid = FLUID;
                var ghost = GHOST;
                var numberofparticles = NumberOfParticles;
                gpu.For(0, NumberOfParticles, i =>
                {
                    if (particletype[i] == fluid)
                    {
                        double viscosityTerm_x = 0.0;
                        double viscosityTerm_y = 0.0;
                        double viscosityTerm_z = 0.0;
                        for (int j = 0; j < numberofparticles; j++)
                        {
                            if ((j == i) || (particletype[j] == ghost)) continue;
                            double xij = position[j * 3] - position[i * 3];
                            double yij = position[j * 3 + 1] - position[i * 3 + 1];
                            double zij = position[j * 3 + 2] - position[i * 3 + 2];
                            double distance2 = (xij * xij) + (yij * yij) + (zij * zij);
                            double distance = Math.Sqrt(distance2);
                            if (distance < re)
                            {
                                double w = re / distance - 1;
                                viscosityTerm_x += (velocity[j * 3] - velocity[i * 3]) * w;
                                viscosityTerm_y += (velocity[j * 3 + 1] - velocity[i * 3 + 1]) * w;
                                viscosityTerm_z += (velocity[j * 3 + 2] - velocity[i * 3 + 2]) * w;
                            }
                        }
                        viscosityTerm_x = viscosityTerm_x * a;
                        viscosityTerm_y = viscosityTerm_y * a;
                        viscosityTerm_z = viscosityTerm_z * a;
                        acceleration[i * 3] += viscosityTerm_x;
                        acceleration[i * 3 + 1] += viscosityTerm_y;
                        acceleration[i * 3 + 2] += viscosityTerm_z;
                    }
                });
                Acceleration = acceleration;
            }
        }


        void moveParticle()
        {
            for (int i = 0; i < NumberOfParticles; i++)
            {
                if (ParticleType[i] == FLUID)
                {
                    Velocity[i * 3] += Acceleration[i * 3] * DT;
                    Velocity[i * 3 + 1] += Acceleration[i * 3 + 1] * DT;
                    Velocity[i * 3 + 2] += Acceleration[i * 3 + 2] * DT;

                    Position[i * 3] += Velocity[i * 3] * DT;
                    Position[i * 3 + 1] += Velocity[i * 3 + 1] * DT;
                    Position[i * 3 + 2] += Velocity[i * 3 + 2] * DT;
                }
                Acceleration[i * 3] = 0.0;
                Acceleration[i * 3 + 1] = 0.0;
                Acceleration[i * 3 + 2] = 0.0;
            }
        }


        void collision()
        {
            int i, j;
            double xij, yij, zij;
            double distance, distance2;
            double forceDT; /* forceDT is the impulse of collision between particles */
            double mi, mj;
            double velocity_ix, velocity_iy, velocity_iz;
            double e = COEFFICIENT_OF_RESTITUTION;
            var VelocityAfterCollision = new double[3 * ARRAY_SIZE]; // 元はstaticになっていた。

            for (i = 0; i < 3 * NumberOfParticles; i++)
            {
                VelocityAfterCollision[i] = Velocity[i];
            }
            for (i = 0; i < NumberOfParticles; i++)
            {
                if (ParticleType[i] == FLUID)
                {
                    mi = FluidDensity;
                    velocity_ix = Velocity[i * 3];
                    velocity_iy = Velocity[i * 3 + 1];
                    velocity_iz = Velocity[i * 3 + 2];
                    for (j = 0; j < NumberOfParticles; j++)
                    {
                        if ((j == i) || (ParticleType[j] == GHOST)) continue;
                        xij = Position[j * 3] - Position[i * 3];
                        yij = Position[j * 3 + 1] - Position[i * 3 + 1];
                        zij = Position[j * 3 + 2] - Position[i * 3 + 2];
                        distance2 = (xij * xij) + (yij * yij) + (zij * zij);
                        if (distance2 < collisionDistance2)
                        {
                            distance = Math.Sqrt(distance2);
                            forceDT = (velocity_ix - Velocity[j * 3]) * (xij / distance)
                                     + (velocity_iy - Velocity[j * 3 + 1]) * (yij / distance)
                                     + (velocity_iz - Velocity[j * 3 + 2]) * (zij / distance);
                            if (forceDT > 0.0)
                            {
                                mj = FluidDensity;
                                forceDT *= (1.0 + e) * mi * mj / (mi + mj);
                                velocity_ix -= (forceDT / mi) * (xij / distance);
                                velocity_iy -= (forceDT / mi) * (yij / distance);
                                velocity_iz -= (forceDT / mi) * (zij / distance);
                                /*
                                if(j>i){ fConsole.WriteLine(stderr,"WARNING: Collision occured between %d and %d particles.\n",i,j); }
                                */
                            }
                        }
                    }
                    VelocityAfterCollision[i * 3] = velocity_ix;
                    VelocityAfterCollision[i * 3 + 1] = velocity_iy;
                    VelocityAfterCollision[i * 3 + 2] = velocity_iz;
                }
            }
            for (i = 0; i < NumberOfParticles; i++)
            {
                if (ParticleType[i] == FLUID)
                {
                    Position[i * 3] += (VelocityAfterCollision[i * 3] - Velocity[i * 3]) * DT;
                    Position[i * 3 + 1] += (VelocityAfterCollision[i * 3 + 1] - Velocity[i * 3 + 1]) * DT;
                    Position[i * 3 + 2] += (VelocityAfterCollision[i * 3 + 2] - Velocity[i * 3 + 2]) * DT;
                    Velocity[i * 3] = VelocityAfterCollision[i * 3];
                    Velocity[i * 3 + 1] = VelocityAfterCollision[i * 3 + 1];
                    Velocity[i * 3 + 2] = VelocityAfterCollision[i * 3 + 2];
                }
            }
        }


        void calPressure()
        {
            calNumberDensity();
            setBoundaryCondition();
            setSourceTerm();
            setMatrix();
            solveSimultaniousEquationsByGaussEliminationMethod();
            removeNegativePressure();
            setMinimumPressure();
        }


        void calNumberDensity()
        {
            if (!usegpu)
            {
                for (int i = 0; i < NumberOfParticles; i++)
                {
                    NumberDensity[i] = 0.0;
                    if (ParticleType[i] == GHOST) continue;
                    for (int j = 0; j < NumberOfParticles; j++)
                    {
                        if ((j == i) || (ParticleType[j] == GHOST)) continue;
                        double xij = Position[j * 3] - Position[i * 3];
                        double yij = Position[j * 3 + 1] - Position[i * 3 + 1];
                        double zij = Position[j * 3 + 2] - Position[i * 3 + 2];
                        double distance2 = (xij * xij) + (yij * yij) + (zij * zij);
                        double distance = Math.Sqrt(distance2);
                        double w = weight(distance, Re_forNumberDensity);
                        NumberDensity[i] += w;
                    }
                }
            }
            else
            {
                var particletype = ParticleType;
                var position = Position;
                var ghost = GHOST;
                var re = Re_forNumberDensity;
                var numberdensity = new double[NumberOfParticles];
                var particlenum = NumberOfParticles;
                gpu.For(0, particlenum, i =>
                {
                    double temp = 0;
                    if (particletype[i] != ghost)
                    {
                        for (int j = 0; j < particlenum; j++)
                        {
                            if ((j == i) || (particletype[j] == ghost)) continue;
                            double xij = position[j * 3] - position[i * 3];
                            double yij = position[j * 3 + 1] - position[i * 3 + 1];
                            double zij = position[j * 3 + 2] - position[i * 3 + 2];
                            double distance2 = (xij * xij) + (yij * yij) + (zij * zij);
                            double distance = Math.Sqrt(distance2);
                            if (distance < re)
                            {
                                temp += re / distance - 1;
                            }
                        }
                    }
                    numberdensity[i] = temp;
                });
                NumberDensity = numberdensity;
            }
        }


        void setBoundaryCondition()
        {
            int i;
            double n0 = N0_forNumberDensity;
            double beta = THRESHOLD_RATIO_OF_NUMBER_DENSITY;

            for (i = 0; i < NumberOfParticles; i++)
            {
                if (ParticleType[i] == GHOST || ParticleType[i] == DUMMY_WALL)
                {
                    BoundaryCondition[i] = GHOST_OR_DUMMY;
                }
                else if (NumberDensity[i] < beta * n0)
                {
                    BoundaryCondition[i] = SURFACE_PARTICLE;
                }
                else
                {
                    BoundaryCondition[i] = INNER_PARTICLE;
                }
            }
        }


        void setSourceTerm()
        {
            int i;
            double n0 = N0_forNumberDensity;
            double gamma = RELAXATION_COEFFICIENT_FOR_PRESSURE;

            for (i = 0; i < NumberOfParticles; i++)
            {
                SourceTerm[i] = 0.0;
                if (ParticleType[i] == GHOST || ParticleType[i] == DUMMY_WALL) continue;
                if (BoundaryCondition[i] == INNER_PARTICLE)
                {
                    SourceTerm[i] = gamma * (1.0 / (DT * DT)) * ((NumberDensity[i] - n0) / n0);
                }
                else if (BoundaryCondition[i] == SURFACE_PARTICLE)
                {
                    SourceTerm[i] = 0.0;
                }
            }
        }


        void setMatrix()
        {
            double xij, yij, zij;
            double distance, distance2;
            double coefficientIJ;
            double n0 = N0_forLaplacian;
            int i, j;
            double a;
            int n = NumberOfParticles;

            for (i = 0; i < NumberOfParticles; i++)
            {
                for (j = 0; j < NumberOfParticles; j++)
                {
                    CoefficientMatrix[i * n + j] = 0.0;
                }
            }

            a = 2.0 * DIM / (n0 * Lambda);
            for (i = 0; i < NumberOfParticles; i++)
            {
                if (BoundaryCondition[i] != INNER_PARTICLE) continue;
                for (j = 0; j < NumberOfParticles; j++)
                {
                    if ((j == i) || (BoundaryCondition[j] == GHOST_OR_DUMMY)) continue;
                    xij = Position[j * 3] - Position[i * 3];
                    yij = Position[j * 3 + 1] - Position[i * 3 + 1];
                    zij = Position[j * 3 + 2] - Position[i * 3 + 2];
                    distance2 = (xij * xij) + (yij * yij) + (zij * zij);
                    distance = Math.Sqrt(distance2);
                    if (distance >= Re_forLaplacian) continue;
                    coefficientIJ = a * weight(distance, Re_forLaplacian) / FluidDensity;
                    CoefficientMatrix[i * n + j] = (-1.0) * coefficientIJ;
                    CoefficientMatrix[i * n + i] += coefficientIJ;
                }
                CoefficientMatrix[i * n + i] += (COMPRESSIBILITY) / (DT * DT);
            }
            exceptionalProcessingForBoundaryCondition();
        }


        void exceptionalProcessingForBoundaryCondition()
        {
            /* If tere is no Dirichlet boundary condition on the fluid, 
               increase the diagonal terms of the matrix for an exception. This allows us to solve the matrix without Dirichlet boundary conditions. */
            checkBoundaryCondition();
            increaseDiagonalTerm();
        }


        void checkBoundaryCondition()
        {
            int i, j, count;
            double xij, yij, zij, distance2;

            for (i = 0; i < NumberOfParticles; i++)
            {
                if (BoundaryCondition[i] == GHOST_OR_DUMMY)
                {
                    FlagForCheckingBoundaryCondition[i] = GHOST_OR_DUMMY;
                }
                else if (BoundaryCondition[i] == SURFACE_PARTICLE)
                {
                    FlagForCheckingBoundaryCondition[i] = DIRICHLET_BOUNDARY_IS_CONNECTED;
                }
                else
                {
                    FlagForCheckingBoundaryCondition[i] = DIRICHLET_BOUNDARY_IS_NOT_CONNECTED;
                }
            }

            do
            {
                count = 0;
                for (i = 0; i < NumberOfParticles; i++)
                {
                    if (FlagForCheckingBoundaryCondition[i] == DIRICHLET_BOUNDARY_IS_CONNECTED)
                    {
                        for (j = 0; j < NumberOfParticles; j++)
                        {
                            if (j == i) continue;
                            if ((ParticleType[j] == GHOST) || (ParticleType[j] == DUMMY_WALL)) continue;
                            if (FlagForCheckingBoundaryCondition[j] == DIRICHLET_BOUNDARY_IS_NOT_CONNECTED)
                            {
                                xij = Position[j * 3] - Position[i * 3];
                                yij = Position[j * 3 + 1] - Position[i * 3 + 1];
                                zij = Position[j * 3 + 2] - Position[i * 3 + 2];
                                distance2 = (xij * xij) + (yij * yij) + (zij * zij);
                                if (distance2 >= Re2_forLaplacian) continue;
                                FlagForCheckingBoundaryCondition[j] = DIRICHLET_BOUNDARY_IS_CONNECTED;
                            }
                        }
                        FlagForCheckingBoundaryCondition[i] = DIRICHLET_BOUNDARY_IS_CHECKED;
                        count++;
                    }
                }
            } while (count != 0); /* This procedure is repeated until the all fluid or wall particles (which have Dirhchlet boundary condition in the particle group) are in the state of "DIRICHLET_BOUNDARY_IS_CHECKED".*/

            for (i = 0; i < NumberOfParticles; i++)
            {
                if (FlagForCheckingBoundaryCondition[i] == DIRICHLET_BOUNDARY_IS_NOT_CONNECTED)
                {
                    Console.WriteLine("WARNING: There is no dirichlet boundary condition for %d-th particle.\n", i); // 元はfprintfでstderrに書き込んでいた
                }
            }
        }


        void increaseDiagonalTerm()
        {
            int i;
            int n = NumberOfParticles;

            for (i = 0; i < n; i++)
            {
                if (FlagForCheckingBoundaryCondition[i] == DIRICHLET_BOUNDARY_IS_NOT_CONNECTED)
                {
                    CoefficientMatrix[i * n + i] = 2.0 * CoefficientMatrix[i * n + i];
                }
            }
        }


        void solveSimultaniousEquationsByGaussEliminationMethod()
        {
            // ほぼすべての時間がこれの計算にかかっている
            int n = NumberOfParticles;
            if (!usegpu)
            {
                for (int i = 0; i < n; i++)
                {
                    Pressure[i] = 0.0;
                }
                for (int i = 0; i < n - 1; i++)
                {
                    if (BoundaryCondition[i] != INNER_PARTICLE) continue;
                    for (int j = i + 1; j < n; j++)
                    {
                        if (BoundaryCondition[j] == GHOST_OR_DUMMY) continue;
                        double c = CoefficientMatrix[j * n + i] / CoefficientMatrix[i * n + i];
                        for (int k = i + 1; k < n; k++)
                        {
                            CoefficientMatrix[j * n + k] -= c * CoefficientMatrix[i * n + k];
                        }
                        SourceTerm[j] -= c * SourceTerm[i];
                    }
                }
                for (int i = n - 1; i >= 0; i--)
                {
                    if (BoundaryCondition[i] != INNER_PARTICLE) continue;
                    double sumOfTerms = 0.0;
                    for (int j = i + 1; j < n; j++)
                    {
                        if (BoundaryCondition[j] == GHOST_OR_DUMMY) continue;
                        sumOfTerms += CoefficientMatrix[i * n + j] * Pressure[j];
                    }
                    Pressure[i] = (SourceTerm[i] - sumOfTerms) / CoefficientMatrix[i * n + i];
                }
            }
            else
            {
                Pressure = ConjugateGradientMethod();
            }
        }


        void removeNegativePressure()
        {
            int i;

            for (i = 0; i < NumberOfParticles; i++)
            {
                if (Pressure[i] < 0.0) Pressure[i] = 0.0;
            }
        }


        void setMinimumPressure()
        {
            if (!usegpu)
            {
                for (int i = 0; i < NumberOfParticles; i++)
                {
                    if (ParticleType[i] == GHOST || ParticleType[i] == DUMMY_WALL) continue;
                    MinimumPressure[i] = Pressure[i];
                    for (int j = 0; j < NumberOfParticles; j++)
                    {
                        if ((j == i) || (ParticleType[j] == GHOST)) continue;
                        if (ParticleType[j] == DUMMY_WALL) continue;
                        double xij = Position[j * 3] - Position[i * 3];
                        double yij = Position[j * 3 + 1] - Position[i * 3 + 1];
                        double zij = Position[j * 3 + 2] - Position[i * 3 + 2];
                        double distance2 = (xij * xij) + (yij * yij) + (zij * zij);
                        if (distance2 >= Re2_forGradient) continue;
                        if (MinimumPressure[i] > Pressure[j])
                        {
                            MinimumPressure[i] = Pressure[j];
                        }
                    }
                }
            }
            else
            {
                var particletype = ParticleType;
                var position = Position;
                var pressure = Pressure;
                var n = NumberOfParticles;
                var minpressure = new double[n];
                var dw = DUMMY_WALL;
                var ghost = GHOST;
                var re = Re2_forGradient;
                gpu.For(0, n, i =>
                {
                    if (!(particletype[i] == ghost || particletype[i] == dw))
                    {
                        double temp = pressure[i];
                        for (int j = 0; j < n; j++)
                        {
                            if ((j == i) || (particletype[j] == ghost)) continue;
                            if (particletype[j] == dw) continue;
                            double xij = position[j * 3] - position[i * 3];
                            double yij = position[j * 3 + 1] - position[i * 3 + 1];
                            double zij = position[j * 3 + 2] - position[i * 3 + 2];
                            double distance2 = (xij * xij) + (yij * yij) + (zij * zij);
                            if (distance2 >= re) continue;
                            if (temp > pressure[j])
                            {
                                temp = pressure[j];
                            }
                        }
                        minpressure[i] = temp;
                    }
                });
                MinimumPressure = minpressure;
            }
        }


        void calPressureGradient()
        {
            double a = DIM / N0_forGradient;
            if (!usegpu)
            {
                for (int i = 0; i < NumberOfParticles; i++)
                {
                    if (ParticleType[i] != FLUID) continue;
                    double gradient_x = 0.0;
                    double gradient_y = 0.0;
                    double gradient_z = 0.0;
                    for (int j = 0; j < NumberOfParticles; j++)
                    {
                        if (j == i) continue;
                        if (ParticleType[j] == GHOST) continue;
                        if (ParticleType[j] == DUMMY_WALL) continue;
                        double xij = Position[j * 3] - Position[i * 3];
                        double yij = Position[j * 3 + 1] - Position[i * 3 + 1];
                        double zij = Position[j * 3 + 2] - Position[i * 3 + 2];
                        double distance2 = (xij * xij) + (yij * yij) + (zij * zij);
                        double distance = Math.Sqrt(distance2);
                        if (distance < Re_forGradient)
                        {
                            double w = weight(distance, Re_forGradient);
                            double pij = (Pressure[j] - MinimumPressure[i]) / distance2;
                            gradient_x += xij * pij * w;
                            gradient_y += yij * pij * w;
                            gradient_z += zij * pij * w;
                        }
                    }
                    gradient_x *= a;
                    gradient_y *= a;
                    gradient_z *= a;
                    Acceleration[i * 3] = (-1.0) * gradient_x / FluidDensity;
                    Acceleration[i * 3 + 1] = (-1.0) * gradient_y / FluidDensity;
                    Acceleration[i * 3 + 2] = (-1.0) * gradient_z / FluidDensity;
                }
            }
            else
            {
                var particletype = ParticleType;
                var position = Position;
                var pressure = Pressure;
                var minpressure = MinimumPressure;
                var fluid = FLUID;
                var ghost = GHOST;
                var dw = DUMMY_WALL;
                int n = NumberOfParticles;
                var re = Re_forGradient;
                var fd = FluidDensity;
                var acceleration = new double[n * 3];
                gpu.For(0, n, i =>
                {
                    if (particletype[i] == fluid)
                    {
                        double gradient_x = 0.0;
                        double gradient_y = 0.0;
                        double gradient_z = 0.0;
                        for (int j = 0; j < n; j++)
                        {
                            if (j == i) continue;
                            if (particletype[j] == ghost) continue;
                            if (particletype[j] == dw) continue;
                            double xij = position[j * 3] - position[i * 3];
                            double yij = position[j * 3 + 1] - position[i * 3 + 1];
                            double zij = position[j * 3 + 2] - position[i * 3 + 2];
                            double distance2 = (xij * xij) + (yij * yij) + (zij * zij);
                            double distance = Math.Sqrt(distance2);
                            if (distance < re)
                            {
                                double w = re / distance - 1;
                                double pij = (pressure[j] - minpressure[i]) / distance2;
                                gradient_x += xij * pij * w;
                                gradient_y += yij * pij * w;
                                gradient_z += zij * pij * w;
                            }
                        }
                        gradient_x *= a;
                        gradient_y *= a;
                        gradient_z *= a;
                        acceleration[i * 3] = (-1.0) * gradient_x / fd;
                        acceleration[i * 3 + 1] = (-1.0) * gradient_y / fd;
                        acceleration[i * 3 + 2] = (-1.0) * gradient_z / fd;
                    }
                });
                Acceleration = acceleration;
            }
        }


        void moveParticleUsingPressureGradient()
        {
            int i;

            for (i = 0; i < NumberOfParticles; i++)
            {
                if (ParticleType[i] == FLUID)
                {
                    Velocity[i * 3] += Acceleration[i * 3] * DT;
                    Velocity[i * 3 + 1] += Acceleration[i * 3 + 1] * DT;
                    Velocity[i * 3 + 2] += Acceleration[i * 3 + 2] * DT;

                    Position[i * 3] += Acceleration[i * 3] * DT * DT;
                    Position[i * 3 + 1] += Acceleration[i * 3 + 1] * DT * DT;
                    Position[i * 3 + 2] += Acceleration[i * 3 + 2] * DT * DT;
                }
                Acceleration[i * 3] = 0.0;
                Acceleration[i * 3 + 1] = 0.0;
                Acceleration[i * 3 + 2] = 0.0;
            }
        }

        private void WriteData()
        {
            positions.Add(Position.ToArray());
            particletypes.Add(ParticleType.ToArray());
        }

        public void ShowData()
        {
            chart.Series[0].Points.Clear();
            chart.Series[1].Points.Clear();
            int frnum = 0;
            for (int i = 0; i < NumberOfParticles; i++)
            {
                if (ParticleType[i] >= 0)
                {
                    int seriesnum = ParticleType[i] < 2 ? 0 : 1;
                    if (seriesnum == 0)
                    {
                        chart.Series[seriesnum].Points.AddXY(Position[i * 3], Position[i * 3 + 1]);
                        chart.Series[seriesnum].Points[frnum].Color = System.Drawing.Color.FromArgb(128, 0, 192, 192);
                        frnum++;
                    }
                    else
                    {
                        chart.Series[seriesnum].Points.AddXY(Position[i * 3], Position[i * 3 + 1]);
                    }
                }
            }
            DoEvent();
        }

        

        private double[] ConjugateGradientMethod()
        {
            int n = NumberOfParticles;
            var weight = CoefficientMatrix;
            var product = SourceTerm;
            var bc = BoundaryCondition;
            var inner = INNER_PARTICLE;
            var ghostordummy = GHOST_OR_DUMMY;
            var befpressure = Pressure;

            var tempijs = new List<int>();
            var tempis = new List<int>();
            var tempjs = new List<int>();
            // Ax = b において jがxに、iがbに対応
            // weightはiが縦、jが横
            // p, restに対応するインデックスはi
            for (int j = 0; j < n; j++)
            {
                if (bc[j] == inner)
                {
                    for (int i = 0; i < n; i++)
                    {
                        if (bc[i] != ghostordummy && weight[i * n + j] != 0)
                        {
                            tempijs.Add(i * n + j);
                            tempis.Add(i);
                            tempjs.Add(j);
                        }
                    }
                }
            }
            var thinkijs = tempijs.ToArray();
            var thinkis = tempis.ToArray();
            var thinkjs = tempjs.ToArray();
            
            var output = befpressure.ToArray();
            double[] wx = Product(weight, output, thinkijs, thinkis, thinkjs);
            var rest = Subtract(product, wx);
            var p = rest.ToArray();
            double threshold = -1;
            int count = 0;
            while (count < 200)
            {
                double[] wp = Product(weight, p, thinkijs, thinkis, thinkjs);
                double a = Product(rest, p) / Product(p, wp);
                output = Add(output, Product(a, p));
                double[] nextrest = Subtract(rest, Product(a, wp));
                double loss = nextrest.Sum(x => x * x);
                if (threshold < 0)
                {
                    threshold = loss / 100;
                }
                else
                {
                    if (loss < threshold)
                    {
                        break;
                    }
                }
                double b = Product(nextrest, nextrest) / Product(rest, rest);
                p = Add(nextrest, Product(b, p));
                Array.Copy(nextrest, rest, rest.Length);
                count++;
            }
            return output;
        }

        



        // ここから一般的な関数
        private double Product(double[] x, double[] y)
        {
            double output = 0;
            for (int i = 0; i < Math.Min(x.Length, y.Length); i++)
            {
                output += x[i] * y[i];
            }
            return output;
        }
        

        private double[] Product(double[] x, double[] y, int[] thinkijs, int[] thinkis, int[] thinkjs)
        {
            var output = new double[y.Length]; // xはn*n, yはn の長さの配列
            for (int ij = 0; ij < thinkijs.Length; ij++)
            {
                output[thinkjs[ij]] += x[thinkijs[ij]] * y[thinkis[ij]];
            }
            return output;
        }


        private double[] Product(double x, double[] y)
        {
            var output = new double[y.Length];
            for (int i = 0; i < y.Length; i++)
            {
                output[i] = x * y[i];
            }
            return output;
        }

        private double[] Subtract(double[] x, double[] y)
        {
            int tate = Math.Min(x.Length, y.Length);
            var output = new double[tate];
            for (int i = 0; i < tate; i++)
            {
                output[i] = x[i] - y[i];
            }
            return output;
        }

        private double[] Add(double[] x, double[] y)
        {
            int tate = Math.Min(x.Length, y.Length);
            var output = new double[tate];
            for (int i = 0; i < tate; i++)
            {
                output[i] = x[i] + y[i];
            }
            return output;
        }

        // ここまで一般的な関数

    }
}

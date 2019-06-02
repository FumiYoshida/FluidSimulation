using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace FluidSimulation
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        private void button2_Click(object sender, EventArgs e)
        {

        }

        private void button1_Click(object sender, EventArgs e)
        {
            var particlenum = 100;
            var particleinfosize = 4;
            var data = new double[particlenum, particleinfosize];
            var rd = new Random();
            for (int i = 0; i < particlenum; i++)
            {
                data[i, 0] = rd.Next(100);
                data[i, 1] = rd.Next(100);
            }

            WriteData(data);


        }

        private void WriteData(double[,] data)
        {
            chart1.Series[0].Points.Clear();
            for (int i = 0; i < data.GetLength(0); i++)
            {
                chart1.Series[0].Points.AddXY(data[i, 0], data[i, 1]);
            }
        }


    }
}

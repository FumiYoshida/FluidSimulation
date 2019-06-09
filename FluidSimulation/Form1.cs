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
        
        ApproxCalc apc;
        private void button1_Click(object sender, EventArgs e)
        {
            apc = new ApproxCalc();
            WriteData(apc);
        }


        private void button2_Click(object sender, EventArgs e)
        {
            label1.Text = "calculating...";
            for (int i = 0; i < 10; i++)
            {
                System.Windows.Forms.Application.DoEvents();
                for (int j = 0; j < 1; j++)
                {
                    apc.Next();
                    WriteData(apc);
                }
            }
            label1.Text = "finished";
        }


        private void WriteData(ApproxCalc apc)
        {
            var data = apc.data[0];
            chart1.Series[0].Points.Clear();
            var pressures = new double[data.GetLength(0)];
            int serial = 0;
            if (apc.pressure != null)
            {
                pressures = apc.pressure;
            }
            for (int i = 0; i < data.GetLength(0); i++)
            {
                if (apc.isdummy[i] | apc.iswall[i])
                {
                    chart1.Series[1].Points.AddXY(data[i, 0], data[i, 1]);
                }
                else
                {
                    chart1.Series[0].Points.AddXY(data[i, 0], data[i, 1]);
                    int cl = (int)(Math.Max(0, Math.Min(255, 128 * (pressures[i] + 0.5))));
                    chart1.Series[0].Points[serial].Color = Color.FromArgb(cl, 0, 0, 0);
                    serial++;
                }
                
            }
        }

        private void chart1_Click(object sender, EventArgs e)
        {

        }

        private void button3_Click(object sender, EventArgs e)
        {
            var mps = new Mps(chart1, System.Windows.Forms.Application.DoEvents, true);
            mps.debuglabel = label1;
            mps.Main();
        }
    }
}

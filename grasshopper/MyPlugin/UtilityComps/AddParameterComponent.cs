using System;
using System.Collections.Generic;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace MyPlugin.UtilityComps
{
    public class AddParameterComponent : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MyComponent2 class.
        /// </summary>
        public AddParameterComponent()
          : base("Add Params", "Params",
              "Description",
              "MyPlugin", "Utilites")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddBooleanParameter("run", "r", "delete the params", GH_ParamAccess.item, false);
            pManager.AddNumberParameter("num", "n", "trying", GH_ParamAccess.list);
            Params.Input[1].Optional = true;
            pManager.AddNumberParameter("v", "n", "trying", GH_ParamAccess.item);
            Params.Input[2].Optional = true;
            pManager.AddBooleanParameter("Reset", "R", "delete the params", GH_ParamAccess.item, false);
            
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddTextParameter("out", "o", "output info", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        /// Using Mahdiyar June 2019 post on https://discourse.mcneel.com/t/c-creating-automatic-sliders/84502
        protected override void SolveInstance(IGH_DataAccess DA)
        { 
            // <Custom additional code>   
            Boolean run = false;
            Boolean reset = false;
            if (!DA.GetData(0, ref run)) return;
            if (!DA.GetData(Params.Input.Count-1, ref reset)) return;
            if (Params.Input[1].SourceCount > 0)
            {
                List<IGH_Param> sources = new List<IGH_Param>(Params.Input[1].Sources);
                //get active GH document with this.OnPingDocument https://www.grasshopper3d.com/forum/topics/get-current-grasshopperdocument-in-c
                this.OnPingDocument().RemoveObjects(sources, false);
            }
            //run should be a button
            //could make run a button on the component itself. https://www.youtube.com/watch?v=NDJK5ljmsF4
            //could make a run button next to each param on the connector component

            if (run)
            {
                //loop through all automated inputs (run and reset are at the first and least indices so skip those)
                for (int i = 1; i < Params.Input.Count - 1; i++)
                {
                    //instantiate  new slider
                    Grasshopper.Kernel.Special.GH_NumberSlider slid = new Grasshopper.Kernel.Special.GH_NumberSlider();
                    slid.CreateAttributes(); //sets up default values, and makes sure your slider doesn't crash rhino

                    //customize slider
                    int inputcount = Params.Input[i].SourceCount;
                    slid.Attributes.Pivot = new System.Drawing.PointF((float)Attributes.DocObject.Attributes.Bounds.Left - slid.Attributes.Bounds.Width - 30, (float)Params.Input[i].Attributes.Bounds.Y + inputcount * 30);
                    slid.Slider.Maximum = 120;
                    slid.Slider.Minimum = 0;
                    slid.Slider.DecimalPlaces = 2;
                    slid.SetSliderValue((decimal)(15));

                    //Until now, the slider is a hypothetical object.
                    // This command makes it 'real' and adds it to the canvas.
                    this.OnPingDocument().AddObject(slid, false);
                    //if there is already a slider there then delete  it with the new slider
                    if (Params.Input[i].SourceCount > 0)
                    {
                        Params.Input[i].RemoveAllSources();
                    }
                    //Connect the new slider to this component
                    Params.Input[i].AddSource(slid);
                }
            }

            //if reset is run then delete all the connected params
            if (reset) {
                //loop through all automated inputs (run and reset are at the first and least indices so skip those)
                for (int i = 1; i < Params.Input.Count-1; i++) {
                    Params.Input[i].RemoveAllSources();
                    //TO DO select all the sources being removed and delete them. currently this only un wires them
                    //either get their ids and remove or select all and delete
                    
                }
            }



        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return null;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("D76FE870-A0CB-42E6-AE9D-AA375DB9B2D0"); }
        }
    }
}
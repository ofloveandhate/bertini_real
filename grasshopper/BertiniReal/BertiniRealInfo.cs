using Grasshopper;
using Grasshopper.Kernel;
using System;
using System.Drawing;

//meta data for plugin

namespace BertiniReal
{
    public class BertiniRealInfo : GH_AssemblyInfo
    {
        public override string Name => "BertiniReal";

        //Return a 24x24 pixel bitmap to represent this GHA library.
        public override Bitmap Icon
        {
            get
            {
                //Return a 24x24 pixel bitmap to represent this GHA library.
                return Properties.Resources.plugin_icon;
            }
        }

    //Return a short string describing the purpose of this GHA library.
    public override string Description => "Plugs, sockets, and manipulations for singular algebraic surfaces computed by bertini_real";

        public override Guid Id => new Guid("cf30f496-469e-4187-9e55-cfa50a0cd566");

        //Return a string identifying you or your company.
        public override string AuthorName => "Caden Joergens and silviana amethyst";

        //Return a string representing your preferred contact details.
        public override string AuthorContact => "https://silviana.org";
    }
}
using Grasshopper.Kernel;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace bertini_real
{
    /* Set some properties for the plugin tab*/
    public class TabProperties : GH_AssemblyPriority
    {
        public override GH_LoadingInstruction PriorityLoad()
        {
            /*Register our plug in icon*/
            var server = Grasshopper.Instances.ComponentServer;
            server.AddCategoryShortName("bertini_real", "br");
            server.AddCategorySymbolName("bertini_real", 'B');
            server.AddCategoryIcon("bertini_real", IconLoader.GetIcon("bertini_real.png"));

            return GH_LoadingInstruction.Proceed;
        }
    }
}

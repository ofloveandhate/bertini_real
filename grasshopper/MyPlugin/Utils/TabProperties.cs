using Grasshopper.Kernel;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MyPlugin.Utils
{
    /* Set some properties for the plugin tab*/
    public class TabProperties : GH_AssemblyPriority
    {
        public override GH_LoadingInstruction PriorityLoad()
        {
            /*Register our plug in icon*/
            var server = Grasshopper.Instances.ComponentServer;
            server.AddCategoryShortName("MyPlugin", "MP");
            server.AddCategorySymbolName("MyPlugin", 'P');
            server.AddCategoryIcon("MyPlugin", Properties.Resources.telephone_icon);

            return GH_LoadingInstruction.Proceed;
        }
    }
}

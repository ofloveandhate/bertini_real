using Grasshopper.Kernel;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BertiniReal.Utils
{
    /* Set some properties for the plugin tab*/
    public class TabProperties : GH_AssemblyPriority
    {
        public override GH_LoadingInstruction PriorityLoad()
        {
            /*Register our plug in icon*/
            var server = Grasshopper.Instances.ComponentServer;
            server.AddCategoryShortName("BertiniReal", "BR");
            server.AddCategorySymbolName("BertiniReal", 'B');
            server.AddCategoryIcon("BertiniReal", Properties.Resources.plugin_icon);

            return GH_LoadingInstruction.Proceed;
        }
    }
}

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
            server.AddCategoryShortName("BertiniReal", "MP");
            server.AddCategorySymbolName("BertiniReal", 'P');
            server.AddCategoryIcon("BertiniReal", Properties.Resources.telephone_icon);

            return GH_LoadingInstruction.Proceed;
        }
    }
}

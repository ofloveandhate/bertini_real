function change_line_width(br_plotter,src,hndl)



prompt = {'Enter new line width:'};
			dlg_title = 'LineWidth';
			num_lines = 1;
			default = {num2str(br_plotter.options.linewidth)};

			output = inputdlg(prompt,dlg_title,num_lines,default);

			if isempty(output)
				return;
			end
			
			new_size = str2num(output{1});
			
			if new_size~=br_plotter.options.linewidth
				
                
                
                switch br_plotter.dimension
                    case {2}
                        f = fieldnames(br_plotter.handles.curves.raw);
                        for ii = 1:length(f)
                            set(br_plotter.handles.curves.raw.(f{ii}),'LineWidth',new_size);
                        end
                        
                        f = fieldnames(br_plotter.handles.curves.refinements);
                        for ii = 1:length(f)
                            set(br_plotter.handles.curves.refinements.(f{ii}),'LineWidth',new_size);
                        end
                        
                    case {1}
                        
                        set(br_plotter.handles.edges,'LineWidth',new_size);
						set(br_plotter.handles.sample_edges,'LineWidth',new_size);
                    otherwise
                        
                end
				
				set(br_plotter.buttons.linewidth,'String',sprintf('LineWidth %i',new_size));
				br_plotter.options.linewidth = new_size;
			end
			
			
			
end
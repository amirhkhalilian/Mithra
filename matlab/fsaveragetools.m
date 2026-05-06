classdef fsaveragetools
    % a class to load freesurfer files and convert T1 coord
    % to fs-average
    properties
        HS % hemi
        fn_fsaverage_pial % file path to fs-average pial surface
        fn_fsaverage_sphere % file path to fs-average sphere surface
        fn_subj_pial % file path to subjcet pial surface
        fn_subj_sphere % file path to subject sphere surface
    end

    methods(Static)
        % fucntions for reading freesurfer pial files
        function [retval] = fread3(fid)
            % read a 3 byte int from file
            % copied from freesurfer

            % Original Author: Bruce Fischl
            % CVS Revision Info:
            %    $Author: nicks $
            %    $Date: 2007/01/10 22:55:09 $
            %    $Revision: 1.2 $
            b1 = fread(fid, 1, 'uchar') ;
            b2 = fread(fid, 1, 'uchar') ;
            b3 = fread(fid, 1, 'uchar') ;
            retval = bitshift(b1, 16) + bitshift(b2,8) + b3 ;
        end

        function [vertex_coords, faces] = read_geometry(fname)
            % 
            % [vertex_coords, faces] = read_surf(fname)
            % reads a the vertex coordinates and face lists from a surface file
            % note that reading the faces from a quad file can take a very long
            % time due to the goofy format that they are stored in. If the faces
            % output variable is not specified, they will not be read so it
            % should execute pretty quickly.

            % original file: read_surf.m in ntools

            TRIANGLE_FILE_MAGIC_NUMBER =  16777214 ;
            QUAD_FILE_MAGIC_NUMBER =  16777215 ;

            fid = fopen(fname, 'rb', 'b') ;
            if (fid < 0)
                str = sprintf('could not open curvature file %s.', fname) ;
                error(str) ;
            end
            magic = fsaveragetools.fread3(fid) ;
            if (magic == QUAD_FILE_MAGIC_NUMBER)
                vnum = fsaveragetools.fread3(fid) ; %number of vertices
                fnum = fsaveragetools.fread3(fid) ; %number of faces
                vertex_coords = fread(fid, vnum*3, 'int16') ./ 100 ;
                if (nargout > 1)
                     for i=1:fnum
                         for n=1:4
                            faces(i,n) = fsaveragetools.fread3(fid);
                        end
                    end
                end
            elseif (magic == TRIANGLE_FILE_MAGIC_NUMBER)
                fgets(fid) ;
                fgets(fid) ;
                vnum = fread(fid, 1, 'int32') ; %number of vertices
                fnum = fread(fid, 1, 'int32') ; %number of faces
                vertex_coords = fread(fid, vnum*3, 'float32') + 1 ;
                %The next line was added by Andrew Dykstra (HST-SHBT), Feb 12, 2009, in
                %order to give 3D coordinates for the vertices
                vertex_coords = reshape(vertex_coords, 3, vnum)';
                faces = fread(fid, fnum*3, 'int32') + 1 ;
                faces = reshape(faces, 3, fnum)' ;
            end
            fclose(fid) ;
        end

        function [curv, fnum] = read_curv(fname)
            % [curv, fnum] = read_curv(fname)
            % reads a binary curvature file into a vector
            % https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
            %
            % Reporting: freesurfer@nmr.mgh.harvard.edu
            %
            % open it as a big-endian file
            fid = fopen(fname, 'rb', 'b') ;
            if (fid < 0)
                 str = sprintf('could not open curvature file %s', fname) ;
                 error(str) ;
            end
            vnum = fsaveragetools.fread3(fid) ;
            NEW_VERSION_MAGIC_NUMBER = 16777215;
            if (vnum == NEW_VERSION_MAGIC_NUMBER)
              vnum = fread(fid, 1, 'int32') ;
              fnum = fread(fid, 1, 'int32') ;
              vals_per_vertex = fread(fid, 1, 'int32') ;
              curv = fread(fid, vnum, 'float') ; 
              fclose(fid) ;
            else
              fnum = fsaveragetools.fread3(fid) ;
              curv = fread(fid, vnum, 'int16') ./ 100 ; 
              fclose(fid) ;
            end
        end

        function [verts, faces] = load_surface(file_name);
            % wrapper function for reading a surface in Matlab
            % if file is in .mat format uses load
            % if file is in {.pial, .sphere, .reg, .inflated} 
            % uses ^read_geometry. 

            [filepath,name,ext] = fileparts(file_name);
            if ismember(ext, {'.pial', '.sphere', '.reg', '.inflated'})
                if nargout==1
                    [verts, ~] = fsaveragetools.read_geometry(file_name);
                    verts = verts -1;
                elseif nargout==2
                    [verts, faces] = fsaveragetools.read_geometry(file_name);
                    verts = verts -1;
                else
                    error('too many output args!')
                end
            elseif strcmp(ext, '.mat')
                if nargout==1
                    verts = load(file_name, 'coords');
                elseif nargout==2
                    brain = load(file_name);
                    faces = brain.faces;
                    verts = brain.coords;
                else
                    error('too many output args!')
                end
            else
                error(sprintf('unsupported file extention: %s', ext))
            end
        end
    end

    methods
        function [obj] = fsaveragetools(varargin)
            p = inputParser;
            addParameter(p, 'HS', 'lh')
            addParameter(p, 'fn_fsaverage_pial', []);
            addParameter(p, 'fn_fsaverage_sphere', []);
            addParameter(p, 'fn_subj_pial', []);
            addParameter(p, 'fn_subj_sphere', []);
            parse(p, varargin{:});
            obj.fn_fsaverage_pial = p.Results.fn_fsaverage_pial;
            obj.fn_fsaverage_sphere = p.Results.fn_fsaverage_sphere;
            obj.fn_subj_pial = p.Results.fn_subj_pial;
            obj.fn_subj_sphere = p.Results.fn_subj_sphere;
            if isempty(obj.fn_fsaverage_pial)
                obj.fn_fsaverage_pial = sprintf('./SampleData/fsaverage/%s.pial',obj.HS);
            end
            if isempty(obj.fn_fsaverage_sphere)
                obj.fn_fsaverage_sphere = sprintf('./SampleData/fsaverage/%s.sphere.reg',obj.HS);
            end
        end

        function [fs_coords, fs_inds] = convert_T1_to_fsaverage(obj, ElecLocs, varargin)
            %converts the given coordinates in subject T1 to fsaverage
            %by first finding the closest vertex on subjuct pial
            %then converting the subject pial to subject sphere
            %then find the closest fsaverage vert on sphere
            %then convert fsaverage sphere to fsaverage pail

            p = inputParser;
            addParameter(p, 'fs_pial_verts', []);
            addParameter(p, 'fs_sphr_verts', []);
            addParameter(p, 'subj_pial_verts', []);
            addParameter(p, 'subj_sphr_verts', []);
            parse(p, varargin{:});
            % pial verts for fs-average
            if isempty(p.Results.fs_pial_verts)
                fs_pial_verts = fsaveragetools.load_surface(obj.fn_fsaverage_pial);
            else
                fs_pial_verts = p.Results.fs_pial_verts;
            end
            % sphere verts for fs-average
            if isempty(p.Results.fs_sphr_verts)
                fs_sphr_verts = fsaveragetools.load_surface(obj.fn_fsaverage_sphere);
            else
                fs_sphr_verts = p.Results.fs_sphr_verts;
            end
            % pial verts for subj
            if isempty(p.Results.subj_pial_verts)
                subj_pial_verts = fsaveragetools.load_surface(obj.fn_subj_pial);
            else
                subj_pial_verts = p.Results.subj_pial_verts;
            end
            % sphere verts for subj
            if isempty(p.Results.subj_sphr_verts)
                subj_sphr_verts = fsaveragetools.load_surface(obj.fn_subj_sphere);
            else
                subj_sphr_verts = p.Results.subj_sphr_verts;
            end

            % init fs_coords
            fs_coords = zeros(size(ElecLocs));
            fs_inds = zeros(size(ElecLocs,1),1);
            % loop over elecs
            for i = 1:size(ElecLocs,1)
                % get closest subj pial vert to ElecLocs
                dist_all = sum((subj_pial_verts-ElecLocs(i,:)).^2,2);
                [dist, subj_ind] = min(dist_all);
                % get elec on subj sphere
                ElecLocSubjSphr = subj_sphr_verts(subj_ind,:);
                % find closest vert on fs-average sphere
                dist_all = sum((fs_sphr_verts-ElecLocSubjSphr).^2,2);
                [dist, fs_ind] = min(dist_all);
                % set the ith coords
                fs_inds(i) = fs_ind;
                fs_coords(i,:) = fs_pial_verts(fs_ind,:);
            end
        end
    end

end


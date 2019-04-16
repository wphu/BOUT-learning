"""
Module to allow BOUT.inp files to be read into python and
manipulated with ease.


Nick Walkden, June 2015
nick.walkden@ccfe.ac.uk
"""


from copy import copy
import os

class BOUTOptions(object):
        """
        Class to store and interact with options from BOUT++

        instantiate with
                myOpts = BOUTOptions()
                myOpts.read_inp('path_to_inp_file')
        or
                myOpts = BOUTOptions('path_to_inp_file')

        To get a list of sections use
                section_list = myOpts.list_sections
        or
                section_list = myOpts.list_sections(verbose=True)
        to print to screen too

        Each section of the input is stored as a dictionary attribute so that,
        if you want all the settings in the section [ddx]

                ddx_opt_dict = myOpts.ddx

        and access individual settings by

                ddx_setting = myOpts.ddx['first']

        Any settings in BOUT.inp without a section are stored in

                root_dict = myOpts.root
        """


        def __init__(self,inp_path=None):
                """
                Initialize with the path to the BOUT.inp file. If None initialize as empty.
                """

                self._sections = ['root']

                for section in self._sections:
                        super(BOUTOptions,self).__setattr__(section,{})

                if inp_path is not None:
                        self.read_inp(inp_path)

        def read_inp(self,inp_path=''):

                try:
                        inpfile = open(os.path.join(inp_path, 'BOUT.inp'),'r')
                except:
                        raise TypeError("ERROR: Could not read file "+\
                                        os.path.join(inp_path, "BOUT.inp"))

                current_section = 'root'
                inplines = inpfile.read().splitlines()
                # Close the file after use
                inpfile.close()
                for line in inplines:
                        #remove white space
                        line = line.replace(" ","")


                        if len(line) > 0 and line[0] is not '#':
                                #Only read lines that are not comments or blank
                                if '[' in line:
                                        #Section header
                                        section = line.split('[')[1].split(']')[0]
                                        current_section = copy(section)
                                        if current_section not in self._sections:
                                                self.add_section(current_section)

                                elif '=' in line:
                                        #option setting
                                        attribute = line.split('=')[0]
                                        value = line.split('=')[1].split('#')[0]
                                        value = value.replace("\n","")
                                        value = value.replace("\t","")
                                        value = value.replace("\r","")
                                        value = value.replace("\"","")
                                        self.__dict__[copy(current_section)][copy(attribute)] = copy(value)
                                else:
                                        pass

        def add_section(self,section):
                self._sections.append(section)
                super(BOUTOptions,self).__setattr__(section,{})

        def remove_section(self,section):
                if section in self._sections:
                        self._sections.pop(self._sections.index(sections))
                        super(BOUTOptions,self).__delattr__(section)
                else:
                        print("WARNING: Section "+section+" not found.\n")

        def list_sections(self,verbose=False):
                if verbose:
                        print("Sections Contained: \n")
                        for section in self._sections:
                                print("\t"+section+"\n")

                return self._sections

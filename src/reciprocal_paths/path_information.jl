################################################################################
#
#   Some functions to print information on the path
#
################################################################################

# Function to print some information on a path
"""
    printInfo(path::Path)
prints information on a `Path` object `path`. 
# Examples
```julia-repl
julia> printInfo(path)
...
```
"""



function printInfo(path::ReciprocalPath)          #{ReciprocalPoint{3}}
   
        # print the complete path
        println("Path overview:")
        # maybe already abort if no points or only one point in path
        if length(path.points) == 0
            println("(no points defined)")
            return nothing
        elseif length(path.points) == 1
            println("  ($(1)) \"$(point(path,1).label)\" at $(point(path,1).point) (only point in path)")
            return nothing
      end
        # alternate between points and segments
        for i in 1:length(path.points)-1
            # print the point
            println("  ($(i)) \"$(point(path,i).label)\" at $(point(path,i).point)")
            # print the outgoing segment

       end
        # print the last point
        println("  ($(length(path.points))) \"$(point(path,length(path.points)).label)\" at $(point(path,length(path.points)).point)")
   
  end


# export the function
export printInfo


"""
    getPathString(path::Path)
compiles a string that represents the given path by chaining together all point names.
# Examples
```julia-repl
julia> getPathString(path)
"X--Gamma--K--M--X"
```
"""


function getPathString(path::ReciprocalPath) #{ReciprocalPoint{3}}  
    # create a new string
    path_string = ""
    # distinguish by point number
    if length(path.points) == 0
        # no points in path
        path_string = "(no points defined)"
    elseif length(path.points) == 1
        # one point in path
        path_string = "$(point(path,1).label)"
    else
        # alternate between points and segments
        for i in 1:length(path.points)-1
            # add the point and a segment
            path_string = "$(path_string)$(point(path,i).label)--"
       end
    
        # add the last point
        path_string = "$(path_string)$(point(path,length(path.points)).label)"
         end
    # return the string
    return path_string

   end
# export the function
export getPathString

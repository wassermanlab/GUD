// on change of select doc 
// hide unselected 
$(function () {
    $("#filterDocs").change(function () {
        var resource = $("#filterDocs").val()
        if (resource === "All") {
            $( ".resourceRow" ).show();
            // unhide all
        } else {
            // hide all except that selected
            $( ".resourceRow" ).hide();
            $( "." + resource ).show();
        }
    });
})


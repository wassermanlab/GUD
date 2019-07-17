
var resources;
$(function () {
    $.getJSON("http://127.0.0.1:5000/json_docs", function (json) {
        resources = json;
    });
})
// on change of select resource 
// build url
// populate param list
$(function () {
    $("#resourceSelect").change(function () {
        $("#paramRows").empty();
        var resource = $("#resourceSelect").val()
        var params = resources[resource]['PARAMS']
        for (param in params) {
            var newRow = document.createElement("tr");
            newRow.id = param + "Row";
            $("#paramRows").append(newRow);
            
            var check = document.createElement("td");
            check.id = param + "Check";
            // <td>
            //     <div class="form-check">
            //         <input class="form-check-input"
            //             type="checkbox" value="" checked="checked" id="rowCheck">
            //     </div>
            // </td>
            if (params[param]['REQUIRED']) {
                
            } else {

            }
            
            var key = document.createElement("td");
            key.id = param + "Key";
            key.innerHTML = param;
             // <td>
            //     <input type="text" class="form-control"
            //         id="keyInput">
            // </td>
            var value = document.createElement("td");
            value.id = param + "Val";
            var input = document.createElement("input");
            input.type = "text";
            input.classList.add("form-control");
            input.id = param + "Input";
            
            var desc = document.createElement("td");
            desc.id = param + "Desc";
            desc.innerHTML = params[param]['DESCRIPTION'];
            
            var req = document.createElement("td");
            req.innerHTML = params[param]['REQUIRED'];
            req.id = param + "Req";
            
            $("#"+newRow.id).append(check);
            $("#"+newRow.id).append(key);
            $("#"+newRow.id).append(value);
            $("#"+value.id).append(input);
            $("#"+newRow.id).append(desc);
            $("#"+newRow.id).append(req);

        }
    });
})

// on selection of parameter 
// build url

// on deselection of parameter 
// build url

// on change of parameter value
// build url

// on send 

// on copy



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
            newRow.id = param;
            $("#paramRows").append(newRow);
            // create cjeclnpx  box
            var check = document.createElement("td");
            check.id = param + "Check";
            var form_div = document.createElement("div");
            form_div.id = param + "DivCheck";
            form_div.classList.add("form-check");
            var check_input = document.createElement("input");
            check_input.classList.add("form-check-input");
            check_input.id = param + "RowCheck";
            check_input.setAttribute("type", "checkbox");
            // create key box
            var key = document.createElement("td");
            key.id = param + "Key";
            key.innerHTML = param;
            // create value box
            var value = document.createElement("td");
            value.id = param + "Val";
            var input = document.createElement("input");
            input.type = "text";
            input.classList.add("form-control");
            input.id = param + "Input";
            // create description box
            var desc = document.createElement("td");
            desc.id = param + "Desc";
            desc.innerHTML = params[param]['DESCRIPTION'];
            // create required  box
            var req = document.createElement("td");
            req.innerHTML = params[param]['REQUIRED'];
            req.id = param + "Req";
            // add elements to dom
            $("#" + newRow.id).append(check);
            $("#" + check.id).append(form_div);
            $("#" + form_div.id).append(check_input);
            if (params[param]['REQUIRED']) {
                $("#" + check_input.id).prop("checked", true);
                $("#" + check_input.id).prop("disabled", true);
            }
            $("#" + newRow.id).append(key);
            $("#" + newRow.id).append(value);
            $("#" + value.id).append(input);
            $("#" + newRow.id).append(desc);
            $("#" + newRow.id).append(req);
        }
        build_url()
    });
})

// builds url based on selections on page
function build_url() {
    url = ""
    // get selected resource
    var resource = $("#resourceSelect").val()
    if (resource == "select") {
        $("#url").val(url);
        return;
    }
    // for each selected row get key value pair
    var parameters = {};
    var params;
    var check;
    var key;
    var val;
    $("#paramRows").children().each(function (index) {
        params = this.id;
        check = $("#" + params + "RowCheck").prop("checked")
        key = $("#" + params + "Key").text()
        val = $("#" + params + "Input").val()
        if (check) {
            parameters[key] = val;
        }
    });
    // combine and set get request row
    parameters = Object.entries(parameters);
    url = resource + "?" + parameters[0][0] + "=" + parameters[0][1];
    for (var i = 1; i < parameters.length; i++) {
        url = url + "&" + parameters[i][0] + "=" + parameters[i][1];
    }
    $("#url").val(url);
}

// on selection of parameter 
// build url

// on deselection of parameter 
// build url

// on change of parameter value
// build url


// on send 

// on copy


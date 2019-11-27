// create table row 
function create_row(param, params) {
    var newRow = document.createElement("tr");
    newRow.id = param;
    $("#paramRows").append(newRow);
    var form_div = document.createElement("div");
    form_div.id = param + "DivCheck";
    form_div.classList.add("form-check");
    // create key box
    var key = document.createElement("td");
    key.id = param + "Key";
    key.innerHTML = param;
    // create value box
    var value = document.createElement("td");
    value.id = param + "Val";
    if (param === "genome" || param === "chrom") {
        var input = document.createElement("select")
        input.classList.add("form-control");
        input.id = param + "Input";
    } else {
        var input = document.createElement("input");
        input.type = "text";
        input.classList.add("form-control");
        input.id = param + "Input";
    }
    
    // create description box
    var desc = document.createElement("td");
    desc.id = param + "Desc";
    desc.innerHTML = params[param]['DESCRIPTION'];
    // add elements to dom
    $("#" + newRow.id).append(key);
    $("#" + newRow.id).append(value);
    $("#" + value.id).append(input);
    $("#" + newRow.id).append(desc);
}

// on change of select resource 
// build url
// populate param list
$(function () {
    $("#resourceSelect").change(function () {
        $("#paramRows").empty();
        var resource = $("#resourceSelect").val()
        if (resource === "select"){
            build_url()
            return
        }
        var r = JSON.parse(resources)
        var params = r[resource]['PARAMS']
        var keys = Object.keys(params)
        if (keys.includes("genome")) {
            create_row("genome", params)
            $('#genomeInput').append($('<option>', {value: "hg38",text: 'hg38'}));
            $('#genomeInput').append($('<option>', {value: "hg19",text: 'hg19'}));
        }
        if (keys.includes("chrom")) {
            create_row("chrom", params)
            let val;
            $('#chromInput').append($('<option>', {value: "",text: ''}));
            for (let i=1; i<23 ; i++) {
                val = i.toString()
                $('#chromInput').append($('<option>', {value: val, text: val}));
            }
            $('#chromInput').append($('<option>', {value: "X",text: "X"}));
            $('#chromInput').append($('<option>', {value: "Y",text: "Y"}));
            $('#chromInput').append($('<option>', {value: "M",text: "M"}));
        }
        if (keys.includes("start")) {
            create_row("start", params)
        }
        if (keys.includes("end")) {
            create_row("end", params)
        }
        if (keys.includes("location")) {
            create_row("location", params)
        }

        for (param in params) {
            if (!["chrom", "start", "end", "location", "genome"].includes(param)) {
                create_row(param, params)
            }
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
    var key;
    var val;
    $("#paramRows").children().each(function (index) {
        params = this.id;
        key = $("#" + params + "Key").text()
        val = $("#" + params + "Input").val()
        if (val != "") {
            parameters[key] = val;
        }
    });
    
    // get genome only and remove from rest of parameters 
    genome = parameters["genome"]
    delete parameters.genome
    
    // combine and set get request row
    url = resource
    parameters = Object.entries(parameters);
    if (parameters.length != 0) {
        url = url + "?" + parameters[0][0] + "=" + parameters[0][1];
        for (var i = 1; i < parameters.length; i++) {
            url = url + "&" + parameters[i][0] + "=" + parameters[i][1];
        }
    }

    url = url.replace("{genome}", genome)

    $("#url").val(url);
}

// on change of parameter value
// build url
$(document).on('change', '.form-control', function () {
    build_url()
});

// on send 
$(function () {
    $("#sendButton").click(function () {
        url = $("#url").val();
        url = address_base + url;
        $(".responseCode").html("Loading ...")
        $.ajax({
            url: url,
            dataType: 'json',
            success: function( data ) {
                $(".responseCode").html(JSON.stringify(data, null, 2));
            },
            error: function( data ) {
                error = data['responseText'].match(/>(.*?)<\//g);
                json = {status: error[0].substr(1, error[0].length-3),
                     statusText: error[1].substr(1, error[1].length-3),
                     error: error[2].substr(1, error[2].length-3)}
                $(".responseCode").html(JSON.stringify(json, null, 2));
            }
          });



        // $.getJSON(url, function (json) {
        //     $(".responseCode").html(JSON.stringify(json, null, 2))
        // });
    });
});

// on copy
$(function () {
    $("#copyButton").click(function () {
        var copyText = document.querySelector("#url");
        copyText.select();
        document.execCommand("copy");
    });
});


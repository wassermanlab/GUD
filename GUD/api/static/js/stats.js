// on click to get stats
function getStats(){
    var form = $('#statsForm').serializeArray();
    var url = address_base + "/stats/" + form[0]["value"] + "/" + form[1]["value"];
    window.location.replace(url);
}
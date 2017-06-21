package require http

proc http_get {url args} {
	set h [http::geturl $url -query [http::formatQuery $args]]
	set data [http::data $h]
	http::cleanup $h
	return $data
}

"""Sequence upload endpoint."""
from rest_framework import status
from rest_framework.parsers import FormParser, MultiPartParser
from rest_framework.response import Response
from rest_framework.views import APIView

from apps.sequences.parsing import MAX_UPLOAD_BYTES, ParseError, parse_sequence_file


class ParseSequenceFileView(APIView):
    """POST /api/sequences/parse/ — multipart upload of a FASTA/GenBank file.

    Returns the parsed record with viewer-ready annotations. Nothing is
    written to the database: the client decides whether to keep it, and
    saves it through /api/projects/.
    """

    parser_classes = (MultiPartParser, FormParser)

    def post(self, request):
        upload = request.FILES.get("file")
        if upload is None:
            return Response(
                {"detail": "No file was uploaded. Send it as the 'file' field."},
                status=status.HTTP_400_BAD_REQUEST,
            )
        if upload.size > MAX_UPLOAD_BYTES:
            return Response(
                {"detail": f"File is larger than {MAX_UPLOAD_BYTES // (1024 * 1024)} MB."},
                status=status.HTTP_413_REQUEST_ENTITY_TOO_LARGE,
            )
        try:
            record = parse_sequence_file(upload.read(), filename=upload.name)
        except ParseError as exc:
            return Response({"detail": str(exc)}, status=status.HTTP_400_BAD_REQUEST)
        return Response(record.to_dict())
